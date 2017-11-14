#include "parse.h"

#include "nemlerror.h"

#include <sstream>
#include <iostream>

namespace neml {

std::unique_ptr<NEMLModel> parse_xml(std::string fname, std::string mname,
                                     int & ier)
{
  // Default
  ier = SUCCESS;

  // Parse the XML file
  xmlpp::DomParser parser;

  // Parse, catching file errors
  try {
    parser.parse_file(fname);
  }
  catch(xmlpp::internal_error) {
    ier = FILE_NOT_FOUND;
    return std::unique_ptr<NEMLModel>(nullptr);
  }
  
  try {
    // Grab the root node
    const auto root = parser.get_document()->get_root_node();

    // Dispatch to the right node
    return  find_and_dispatch(root, "material", "name", mname, 
                             &make_from_node, ier);
  }
  catch(...) {
    ier = UNKNOWN_ERROR;
    return std::unique_ptr<NEMLModel>(nullptr);
  }

}

std::unique_ptr<NEMLModel> make_from_node(const xmlpp::Element * node, int & ier)
{
  return dispatch_attribute_unique<NEMLModel>(node, "type",
                                          {"smallstrain","kmregion"},
                                          {&process_smallstrain,&process_kmregion},
                                          ier);
}

std::unique_ptr<NEMLModel> process_kmregion(const xmlpp::Element * node, int & ier)
{
  // Basically need:
  //  1) A list of n smallstrain models
  //  2) A list of n-1 cutoffs
  //  3) An elastic model
  //  4) Boltzmann constant, in your units
  //  5) Burgers vector, in your units
  //  6) Reference strain rate, in your units
 
  // 1)
  std::vector<std::shared_ptr<NEMLModel_sd>> models;
  const xmlpp::Element * models_node;
  if (not one_child(node, "models", models_node, ier)) {
    return std::unique_ptr<NEMLModel>(nullptr);
  }
  auto model_nodes = models_node->get_children("model");
  for (auto it = model_nodes.begin(); it != model_nodes.end(); ++it) {
    auto modeli = process_smallstrain(dynamic_cast<const xmlpp::Element*>(*it), ier);
    std::shared_ptr<NEMLModel> modeli_sh = std::move(modeli);
    models.push_back(std::dynamic_pointer_cast<NEMLModel_sd>(modeli_sh));
  }

  // 2)
  std::vector<double> gs = vector_constant(node, "gcut", ier);

  // Check for lengths
  if (models.size() != (gs.size() + 1)) {
    ier = INCOMPATIBLE_KM;
    return std::unique_ptr<NEMLModel>(nullptr);
  }

  // 3)
  const xmlpp::Element * elastic_node;
  if (not one_child(node, "elastic", elastic_node, ier)) {
    return std::unique_ptr<NEMLModel>(nullptr);
  }
  std::shared_ptr<LinearElasticModel> emodel = process_linearelastic(elastic_node, ier);
  
  // 4)
  double k = scalar_constant(node, "kboltz", ier);

  // 5) 
  double b = scalar_constant(node, "burgers", ier);

  // 6)
  double eps0 = scalar_constant(node, "refrate", ier);

  // CTE, if provided
  std::shared_ptr<Interpolate> alpha;
  auto a_nodes = node->get_children("alpha");
  if (a_nodes.size() > 0) {
    alpha = process_alpha(dynamic_cast<const xmlpp::Element*>(a_nodes.front()), ier);
  }
  else {
    alpha = std::shared_ptr<Interpolate>(new ConstantInterpolate(0.0));
  }

  // Return final model
  return std::unique_ptr<NEMLModel>(
      new KMRegimeModel(models, gs, emodel, k, b, eps0, alpha));
}


std::unique_ptr<NEMLModel> process_smallstrain(const xmlpp::Element * node, int & ier)
{
  // Need a elastic node, a plastic node, and optionally creep and thermal
  // expansion nodes
  // Elastic
  std::shared_ptr<LinearElasticModel> emodel;
  const xmlpp::Element * elastic_node;
  if (not one_child(node, "elastic", elastic_node, ier)) {
    return std::unique_ptr<NEMLModel>(nullptr);
  }
  emodel = process_linearelastic(elastic_node, ier);
  if (ier != SUCCESS) return std::unique_ptr<NEMLModel>(nullptr);

  // Optional alpha
  std::shared_ptr<Interpolate> alpha;
  auto a_nodes = node->get_children("alpha");
  if (a_nodes.size() > 0) {
    alpha = process_alpha(dynamic_cast<const xmlpp::Element*>(a_nodes.front()), ier);
  }
  else {
    alpha = std::shared_ptr<Interpolate>(new ConstantInterpolate(0.0));
  }

  // Optional creep node
  bool found_creep = false;
  auto c_nodes = node->get_children("creep");
  std::shared_ptr<CreepModel> cmodel = nullptr;
  if (c_nodes.size() > 0) {
    cmodel = process_creep(
        dynamic_cast<const xmlpp::Element*>(c_nodes.front()), ier);
    found_creep = true;
  }
  
  // Logic here because we treat viscoplasticity differently
  const xmlpp::Element * plastic_node;
  if (not one_child(node, "plastic", plastic_node, ier)) {
    return std::unique_ptr<NEMLModel>(nullptr);
  }
  
  // Select then on independent/dependent
  std::string ft;
  if (not one_attribute(plastic_node, "type", ft, ier)) {
    return std::unique_ptr<NEMLModel>(nullptr);
  }

  std::unique_ptr<NEMLModel_sd> model;
  if (ft == "independent") {
    std::shared_ptr<RateIndependentFlowRule> fr =  process_independent(
        plastic_node, ier);

    model =  std::unique_ptr<SmallStrainRateIndependentPlasticity> (
      new SmallStrainRateIndependentPlasticity(emodel, fr, alpha));
  }
  else if (ft == "dependent") {
    std::shared_ptr<ViscoPlasticFlowRule> fr = process_dependent(plastic_node, ier);
    std::shared_ptr<TVPFlowRule> gfr = std::make_shared<TVPFlowRule>(emodel, fr);
    model = std::unique_ptr<GeneralIntegrator>(new GeneralIntegrator(gfr, alpha));
  }
  else if (ft == "perfect") {
    // Surface
    std::shared_ptr<YieldSurface> ys = dispatch_node(plastic_node, "surface",
                                                     &process_surface, ier);
    // Yield stress
    std::shared_ptr<Interpolate> yss = scalar_param(plastic_node, "yield", ier);
    model = std::unique_ptr<SmallStrainPerfectPlasticity>(
        new SmallStrainPerfectPlasticity(emodel, ys, yss, alpha));
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::unique_ptr<NEMLModel>(nullptr);
  }

  if (found_creep) {
    return std::unique_ptr<NEMLModel>(
        new SmallStrainCreepPlasticity(std::move(model), cmodel, alpha));
  }

  return std::unique_ptr<NEMLModel>(std::move(model));
}

std::shared_ptr<CreepModel> process_creep(const xmlpp::Element * node, int & ier)
{
  return dispatch_attribute<CreepModel>(node, "type", {"j2"},
                                        {&process_j2creep},
                                        ier);
}

std::shared_ptr<CreepModel> process_j2creep(const xmlpp::Element * node, int & ier)
{
  // Need a scalarmodel node
  std::shared_ptr<ScalarCreepRule> scr = dispatch_node(node, "scalarmodel", 
                                                       &process_scalarmodel, ier);

  return std::shared_ptr<J2CreepModel>(new J2CreepModel(scr));
}

std::shared_ptr<ScalarCreepRule> process_scalarmodel(const xmlpp::Element * node, int & ier)
{
  // dispatch on type
  return dispatch_attribute<ScalarCreepRule>(node, "type", {"powerlaw", "nortonbailey", "regionkm"},
                                             {&process_powerlaw_creep, &process_nb_creep, &process_regionkm_creep},
                                             ier);
}

std::shared_ptr<ScalarCreepRule> process_powerlaw_creep(const xmlpp::Element * node, int & ier)
{
  // Properties are A and n
  return std::shared_ptr<ScalarCreepRule>(new PowerLawCreep(
          scalar_param(node, "A", ier),
          scalar_param(node, "n", ier)
          ));
}

std::shared_ptr<ScalarCreepRule> process_regionkm_creep(const xmlpp::Element * node, int & ier)
{
  // Must have an "elastic" node
  const xmlpp::Element * elastic_node;
  if (not one_child(node, "elastic", elastic_node, ier)) {
    return std::unique_ptr<ScalarCreepRule>(nullptr);
  }

  return std::shared_ptr<ScalarCreepRule>(new RegionKMCreep(
          vector_constant(node, "cuts", ier),
          vector_constant(node, "As", ier),
          vector_constant(node, "Bs", ier),
          scalar_constant(node, "kboltz", ier),
          scalar_constant(node, "b", ier),
          scalar_constant(node, "eps0", ier),
          process_linearelastic(elastic_node, ier)
          ));
}

std::shared_ptr<ScalarCreepRule> process_nb_creep(const xmlpp::Element * node, int & ier)
{
  // Properties are A, m, and n
  return std::shared_ptr<ScalarCreepRule>(new NortonBaileyCreep(
          scalar_param(node, "A", ier),
          scalar_param(node, "m", ier),
          scalar_param(node, "n", ier)
          ));
}

std::shared_ptr<Interpolate> process_alpha(
    const xmlpp::Element * node, int & ier)
{
  // The parameter has name value
  return scalar_param(node, "value", ier);
}

std::shared_ptr<LinearElasticModel> process_linearelastic(
    const xmlpp::Element * node, int & ier)
{
  // Switch on the type of model
  return dispatch_attribute<LinearElasticModel>(node, "type", 
                                                {"isotropic"},
                                                {&process_isotropiclinearelastic},
                                                ier);
}

std::shared_ptr<LinearElasticModel> process_isotropiclinearelastic(
    const xmlpp::Element * node, int & ier)
{
  // Need a shear and a bulk modulus OR a Young's and poisson's
  if((node->get_children("shear").size() > 0) and 
     (node->get_children("bulk").size() > 0)) {
    // Shear
    std::shared_ptr<ShearModulus> sm = dispatch_node(node, "shear", 
                                                     &process_shearmodulus, ier);

    // Bulk
    std::shared_ptr<BulkModulus> bm = dispatch_node(node, "bulk", 
                                                     &process_bulkmodulus, ier);
    
    return std::shared_ptr<LinearElasticModel>(
        new IsotropicLinearElasticModel(sm,bm));
  }
  else if((node->get_children("youngs").size() > 0) and 
     (node->get_children("poissons").size() > 0)) {
    // Young's modulus
    std::shared_ptr<YoungsModulus> em = dispatch_node(node, "youngs", 
                                                     &process_youngsmodulus, ier);

    // Poisson's ratio
    std::shared_ptr<PoissonsRatio> vm = dispatch_node(node, "poissons", 
                                                     &process_poissonsratio, ier);
    
    return std::shared_ptr<LinearElasticModel>(
        new IsotropicLinearElasticModel(em,vm));
  }
  else {
    ier = NODE_NOT_FOUND;
    std::cerr << "Expected to find either {bulk, shear} or {youngs,poissons} " <<
        "for isotropic linear elastic model near line " << node->get_line() <<
        std::endl;
    return std::shared_ptr<LinearElasticModel>(nullptr);
  }
}

std::shared_ptr<ShearModulus> process_shearmodulus(
    const xmlpp::Element * node, int & ier)
{
  // Just the scalar parameter
  return std::shared_ptr<ShearModulus>(new ShearModulus(scalar_param(node, "modulus", ier)));
}

std::shared_ptr<BulkModulus> process_bulkmodulus(
    const xmlpp::Element * node, int & ier)
{
  return std::shared_ptr<BulkModulus>(new BulkModulus(scalar_param(node, "modulus", ier)));
}

std::shared_ptr<YoungsModulus> process_youngsmodulus(
    const xmlpp::Element * node, int & ier)
{
  // Just the scalar parameter
  return std::shared_ptr<YoungsModulus>(new YoungsModulus(scalar_param(node, "modulus", ier)));
}

std::shared_ptr<PoissonsRatio> process_poissonsratio(
    const xmlpp::Element * node, int & ier)
{
  return std::shared_ptr<PoissonsRatio>(new PoissonsRatio(scalar_param(node, "modulus", ier)));
}

std::shared_ptr<RateIndependentFlowRule> process_independent(
    const xmlpp::Element * node, int & ier)
{
  // Should have a "rule" node of differing type
  return dispatch_node(node, "rule", &process_rirule, ier);
}

std::shared_ptr<RateIndependentFlowRule> process_rirule(
    const xmlpp::Element * node, int & ier)
{
  // Types are "associative" or "nonassociative_hardening"
  return dispatch_attribute<RateIndependentFlowRule>(node, "type",
                                          {"associative", "nonassociative_hardening"},
                                          {&process_riassociative, &process_rinonassociative_hardening},
                                          ier);
}

std::shared_ptr<RateIndependentFlowRule> process_riassociative(
    const xmlpp::Element * node, int & ier)
{
  // Need a surface and a hardening model
  // Surface
  std::shared_ptr<YieldSurface> ys = dispatch_node(node, "surface",
                                                   &process_surface, ier);
  // Hardening model
  std::shared_ptr<HardeningRule> hr = dispatch_node(node, "hardening", 
                                                    &process_hardening, ier);

  return std::make_shared<RateIndependentAssociativeFlow>(ys, hr);
}

std::shared_ptr<YieldSurface> process_surface(
    const xmlpp::Element * node, int & ier)
{
  // Switch on type: isoj2, isokinj2
  return dispatch_attribute<YieldSurface>(node, "type",
                                          {"isoj2", "isokinj2",
                                          "isoj2i1", "isokinj2i1"},
                                          {&process_isoj2, &process_isokinj2,
                                          &process_isoj2i1, &process_isokinj2i1},
                                          ier);
}

std::shared_ptr<YieldSurface> process_isoj2(const xmlpp::Element * node,
                                            int & ier)
{
  // No parameters!
  return std::make_shared<IsoJ2>();
}

std::shared_ptr<YieldSurface> process_isokinj2(const xmlpp::Element * node,
                                               int & ier)
{
  // No parameters!
  return std::make_shared<IsoKinJ2>();
}

std::shared_ptr<YieldSurface> process_isoj2i1(const xmlpp::Element * node,
                                            int & ier)
{
  auto h = scalar_param(node, "h", ier);
  auto l = scalar_param(node, "l", ier);
  return std::make_shared<IsoJ2I1>(h,l);
}

std::shared_ptr<YieldSurface> process_isokinj2i1(const xmlpp::Element * node,
                                               int & ier)
{
  auto h = scalar_param(node, "h", ier);
  auto l = scalar_param(node, "l", ier);
  return std::make_shared<IsoKinJ2I1>(h,l);
}

std::shared_ptr<HardeningRule> process_hardening(
    const xmlpp::Element * node, int & ier)
{
  // Three options: isotropic, kinematic, or combined
  // combined will actually recurse to this function, the other two split off
  return dispatch_attribute<HardeningRule>(node, "type",
                                          {"isotropic", "kinematic", "combined"},
                                          {&process_isotropic, &process_kinematic, &process_combined},
                                          ier);
}

std::shared_ptr<HardeningRule> process_isotropic(
    const xmlpp::Element * node, int & ier)
{ 
  // Should have an "isotropic" tag with a type
  return dispatch_node(node, "isotropic", &process_isotropictag, ier);
}

std::shared_ptr<HardeningRule> process_isotropictag(
    const xmlpp::Element * node, int & ier)
{
  return dispatch_attribute<HardeningRule>(node, "type",
                                          {"linear", "interpolated", "voce"},
                                          {&process_linearisotropic, &process_interpolatedisotropic, &process_voceisotropic},
                                          ier);
}

std::shared_ptr<HardeningRule> process_linearisotropic(
    const xmlpp::Element * node, int & ier)
{
  // Two parameters: "yield" and "harden"
  std::shared_ptr<Interpolate> yield = scalar_param(node, "yield", ier);

  std::shared_ptr<Interpolate> harden = scalar_param(node, "harden", ier);

  return std::make_shared<LinearIsotropicHardeningRule>(yield, harden);
}

std::shared_ptr<HardeningRule> process_interpolatedisotropic(
    const xmlpp::Element * node, int & ier)
{
  // Just the flow curve
  std::shared_ptr<Interpolate> flow = scalar_param(node, "flow", ier);

  return std::make_shared<InterpolatedIsotropicHardeningRule>(flow);
}

std::shared_ptr<HardeningRule> process_voceisotropic(
    const xmlpp::Element * node, int & ier)
{
  // Three parameters: "yield", "r", and "d"
  std::shared_ptr<Interpolate> yield = scalar_param(node, "yield", ier);

  std::shared_ptr<Interpolate> r = scalar_param(node, "r", ier);

  std::shared_ptr<Interpolate> d = scalar_param(node, "d", ier);

  return std::make_shared<VoceIsotropicHardeningRule>(yield, r, d);
}


std::shared_ptr<HardeningRule> process_kinematic(
    const xmlpp::Element * node, int & ier)
{
  return dispatch_node(node, "kinematic", &process_kinematictag, ier);
}

std::shared_ptr<HardeningRule> process_kinematictag(
    const xmlpp::Element * node, int & ier)
{
  return dispatch_attribute<HardeningRule>(node, "type",
                                          {"linear"},
                                          {&process_linearkinematic},
                                          ier);
}

std::shared_ptr<HardeningRule> process_linearkinematic(
    const xmlpp::Element * node, int & ier)
{
  // One parameter: "harden"
  std::shared_ptr<Interpolate> harden = scalar_param(node, "harden", ier);

  return std::make_shared<LinearKinematicHardeningRule>(harden);
}

std::shared_ptr<HardeningRule> process_combined(const xmlpp::Element * node,
                                                int & ier)
{
  // Recurse to pick up "isotropic" and "kinematic"
  //
  // Casts should be safe because I dispatch to a function that will result in
  // an error if the types don't match

  // isotropic
  std::shared_ptr<HardeningRule> ir = dispatch_node(node, "isotropic",
                                                    &process_isotropictag, ier);
  std::shared_ptr<IsotropicHardeningRule> cir = 
      std::dynamic_pointer_cast<IsotropicHardeningRule>(ir);

  // kinematic
  std::shared_ptr<HardeningRule> kr = dispatch_node(node, "kinematic", 
                                                    &process_kinematictag, ier);
  std::shared_ptr<KinematicHardeningRule> ckr = 
      std::dynamic_pointer_cast<KinematicHardeningRule>(kr);

  return std::make_shared<CombinedHardeningRule>(cir, ckr);
}

std::shared_ptr<RateIndependentFlowRule> process_rinonassociative_hardening(
    const xmlpp::Element * node, int & ier)
{
  // Need a surface and a (nonassociative) hardening rule
  // Surface
  std::shared_ptr<YieldSurface> ys = dispatch_node(node, "surface",
                                                   &process_surface, ier);

  // Hardening
  std::shared_ptr<NonAssociativeHardening> nsr = dispatch_node(
      node, "hardening", &process_nonass_hardening, ier);

  return std::make_shared<RateIndependentNonAssociativeHardening>(
      ys, nsr);
}

std::shared_ptr<NonAssociativeHardening> process_nonass_hardening(
    const xmlpp::Element * node, int & ier)
{
  // Only Chaboche as an option for now
  return dispatch_attribute<NonAssociativeHardening>(
      node, "type", {"chaboche"}, {&process_nonass_chaboche}, ier);
}

std::shared_ptr<NonAssociativeHardening> process_nonass_chaboche(
    const xmlpp::Element * node, int & ier)
{
  // Need three things: isotropic rule, c vector, and a vector of gamma models
  //
  // isotropic rule
  std::shared_ptr<HardeningRule> ir = dispatch_node(node, "isotropic",
                                                    &process_isotropictag, ier);
  std::shared_ptr<IsotropicHardeningRule> cir = 
      std::dynamic_pointer_cast<IsotropicHardeningRule>(ir);

  // c vector
  std::vector<std::shared_ptr<Interpolate>> c = vector_param(node, "c", ier);
  
  // vector of gamma models
  std::vector< std::shared_ptr<GammaModel> > gvec = 
      dispatch_vector_models<GammaModel>(
          node, "gammamodels", "gammamodel", &process_gammamodel, ier);

  auto matches = node->get_children("A");  
  if (matches.size() > 1) {
    std::vector<std::shared_ptr<Interpolate>> A = vector_param(node, "A", ier);
    std::vector<std::shared_ptr<Interpolate>> a = vector_param(node, "a", ier);
    // Create our object
    return std::make_shared<Chaboche>(cir, c, gvec, A, a);
  }
  else {  
    // Create our object
    return std::make_shared<Chaboche>(cir, c, gvec);
  }
}

std::shared_ptr<GammaModel> process_gammamodel(
    const xmlpp::Element * node, int & ier)
{
  // Two types: constant and saturating
  return dispatch_attribute<GammaModel>(
      node, "type", {"constant", "saturating"},
      {&process_constant_gamma, &process_saturating_gamma},
      ier);
}

std::shared_ptr<GammaModel> process_constant_gamma(
    const xmlpp::Element * node, int & ier)
{
  // Single scalar parameter g
  std::shared_ptr<Interpolate> g = scalar_param(node, "g", ier);

  return std::make_shared<ConstantGamma>(g);
}

std::shared_ptr<GammaModel> process_saturating_gamma(
    const xmlpp::Element * node, int & ier)
{
  // Three scalar parameters: gs, g0, beta
  std::shared_ptr<Interpolate> gs = scalar_param(node, "gs", ier);
  std::shared_ptr<Interpolate> g0 = scalar_param(node, "g0", ier);
  std::shared_ptr<Interpolate> beta = scalar_param(node, "beta", ier);

  return std::make_shared<SatGamma>(gs, g0, beta);
}

std::shared_ptr<ViscoPlasticFlowRule> process_dependent(
    const xmlpp::Element * node, int & ier)
{
  // Should have a "rule" node of differing type
  return dispatch_node(node, "rule", &process_rdrule, ier);
}

std::shared_ptr<ViscoPlasticFlowRule> process_rdrule(
    const xmlpp::Element * node, int & ier)
{
  // Types are "associative", "chaboche", or "yaguchigr91"
  return dispatch_attribute<ViscoPlasticFlowRule>(
      node, "type", {"associative", "chaboche", "yaguchigr91"},
      {&process_rd_associative, &process_rd_chaboche, &process_rd_yaguchigr91},
      ier);
}

std::shared_ptr<ViscoPlasticFlowRule> process_rd_yaguchigr91(
    const xmlpp::Element * node, int & ier)
{
  // Nothing but the model, hardwired coefficients
  return std::make_shared<YaguchiGr91FlowRule>();
}

std::shared_ptr<ViscoPlasticFlowRule> process_rd_chaboche(
    const xmlpp::Element * node, int & ier)
{
  // This is nearly the same as the rate-independent version except we
  // also need a model for the fluidity and the rate n
  // Surface
  std::shared_ptr<YieldSurface> ys = dispatch_node(node, "surface",
                                                   &process_surface, ier);

  // Hardening
  std::shared_ptr<NonAssociativeHardening> nsr = dispatch_node(
      node, "hardening", &process_nonass_hardening, ier);

  // Fluidity
  std::shared_ptr<FluidityModel> fm = dispatch_node(node, "fluidity", 
                                                    &process_fluidity, ier);

  // Rate n
  std::shared_ptr<Interpolate> n = scalar_param(node, "n", ier);

  return std::make_shared<ChabocheFlowRule>(ys, nsr, fm, n);

}

std::shared_ptr<FluidityModel> process_fluidity(
    const xmlpp::Element * node, int & ier)
{
  // Right now only a constant fluidity option
  return dispatch_attribute<FluidityModel>(
      node, "type", {"constant", "saturating"}, 
      {&process_constant_fluidity, &process_saturating_fluidity}, ier);
}

std::shared_ptr<FluidityModel> process_constant_fluidity(
    const xmlpp::Element * node, int & ier)
{
  // One property, eta
  std::shared_ptr<Interpolate> eta = scalar_param(node, "eta", ier);

  return std::make_shared<ConstantFluidity>(eta);
}

std::shared_ptr<FluidityModel> process_saturating_fluidity(
    const xmlpp::Element * node, int & ier)
{
  std::shared_ptr<Interpolate> K0 = scalar_param(node, "K0", ier);
  std::shared_ptr<Interpolate> A = scalar_param(node, "A", ier);
  std::shared_ptr<Interpolate> b = scalar_param(node, "b", ier);

  return std::make_shared<SaturatingFluidity>(K0, A, b);
}

std::shared_ptr<ViscoPlasticFlowRule> process_rd_associative(
    const xmlpp::Element * node, int & ier)
{
  // Almost the same as associative rate independent, but we also need the
  // plastic multiplier function and the fluidity
  // Surface
  std::shared_ptr<YieldSurface> ys = dispatch_node(node, "surface",
                                                   &process_surface, ier);
  // Hardening model
  std::shared_ptr<HardeningRule> hr = dispatch_node(node, "hardening", 
                                                    &process_hardening, ier);

  // gmodel
  std::shared_ptr<GFlow> gm = dispatch_node(node, "gmodel", &process_gmodel,
                                             ier);

  // eta
  std::shared_ptr<Interpolate> eta = scalar_param(node, "eta", ier);

  return std::make_shared<PerzynaFlowRule>(ys, hr, gm, eta);
}

std::shared_ptr<GFlow> process_gmodel(const xmlpp::Element * node, int & ier)
{
  // Just option for power law
  return dispatch_attribute<GFlow>(
      node, "type", {"power_law"}, {&process_gmodel_power_law}, ier);
}

std::shared_ptr<GFlow> process_gmodel_power_law(const xmlpp::Element * node,
                                               int & ier)
{
  // Just the parameter n
  std::shared_ptr<Interpolate> n = scalar_param(node, "n", ier);

  return std::make_shared<GPowerLaw>(n);
}



// Helpers
bool one_child(const xmlpp::Node * node, std::string name,
               const xmlpp::Element * & child, int & ier)
{
  auto matches = node->get_children(name);

  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    std::cerr << "Multiple nodes with name " << name << " found near line " 
        << node->get_line() << std::endl;
    return false;
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
    std::cerr << "Node with name " << name << " not found near line " 
        << node->get_line() << std::endl;
    return false;
  }
  else {
    child = dynamic_cast<const xmlpp::Element*>(matches.front());
    return true;
  }

}

bool one_attribute(const xmlpp::Element * node, std::string name,
                   std::string & value, int & ier)
{
  value = node->get_attribute_value(name);

  if (value == "") {
    ier = ATTRIBUTE_NOT_FOUND;
    std::cerr << "Node " << node->get_name() << " does not have attribute " 
        << name << " near line " << node->get_line() << std::endl;
    return false;
  }
  else {
    return true;
  }
}

std::shared_ptr<Interpolate> scalar_param(const xmlpp::Node * node, 
                                          std::string name, int & ier)
{
  auto matches = node->get_children(name);
  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    std::cerr << "Multiple nodes of type " << name << " found near line " 
        << node->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate());
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
    std::cerr << "Required node type " << name << " not found near line " 
        << node->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate());
  }
  auto child = dynamic_cast<const xmlpp::Element*>(matches.front());
  
  // Two cases: raw text implies constant, interpolate node means use an
  // interpolate

  auto matches_inter = child->get_children("interpolate");

  if (matches_inter.size() != 0) {
    return interpolate_node(matches_inter.front(), ier);
  }
  else if (child->has_child_text()) {
    std::string sval = child->FIRST_TEXT_FN()->get_content();
    if (sval == "") {
      ier = BAD_TEXT;
      std::cerr << "Node " << name << " near line " 
          << node->get_line() << " does not appear to store data" << std::endl;
      return std::shared_ptr<Interpolate>(new InvalidInterpolate());
    }
    return std::shared_ptr<Interpolate>(new ConstantInterpolate(std::stod(sval)));
  }
  else {
    ier = BAD_TEXT;
    std::cerr << "Unable to process scalar parameter input near line " <<
        node->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate());
  }

}

double scalar_constant(const xmlpp::Node * node, std::string name, int & ier)
{
  auto matches = node->get_children(name);
  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    std::cerr << "Multiple nodes of type " << name << " found near line " 
        << node->get_line() << std::endl;
    return 0.0;
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
    std::cerr << "Required node type " << name << " not found near line " 
        << node->get_line() << std::endl;
    return 0.0;
  }
  auto child = dynamic_cast<const xmlpp::Element*>(matches.front());

  if (child->has_child_text()) {
    std::string sval = child->FIRST_TEXT_FN()->get_content();
    if (sval == "") {
      ier = BAD_TEXT;
      std::cerr << "Node " << name << " near line " 
          << node->get_line() << " does not appear to store data" << std::endl;
      return 0.0;
    }
    return std::stod(sval);
  }
  else {
    ier = BAD_TEXT;
    std::cerr << "Unable to process scalar parameter input near line " <<
        node->get_line() << std::endl;
    return 0.0;
  }

}

std::shared_ptr<Interpolate> interpolate_node(const xmlpp::Node * node,
                                              int & ier)
{
  if (node->get_name() != "interpolate") {
    ier = NODE_NOT_FOUND;
    std::cerr << "Interpolation node not found near line" << node->get_line()
        << std::endl;
  }

  auto child = dynamic_cast<const xmlpp::Element*>(node);

  return dispatch_attribute<Interpolate>(child, "type", 
                                         {"constant", "polynomial", "piecewise", "mts", "logpiecewise"},
                                         {&process_constant, &process_polynomial, &process_piecewise, &process_mts, &process_logpiecewise},
                                         ier);
}

std::shared_ptr<Interpolate> process_constant(const xmlpp::Element * child,
                                              int & ier)
{
  // The text is the value
  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    std::cerr << "Constant interpolate node has no text near line " 
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate());
  }
  std::string sval = child->FIRST_TEXT_FN()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    std::cerr << "Constant interpolate node has invalid text near line "
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate()); 
  }
  return std::shared_ptr<Interpolate>(new ConstantInterpolate(std::stod(sval)));
}

std::shared_ptr<Interpolate> process_polynomial(const xmlpp::Element * child,
                                              int & ier)
{
  // The text is a list of coefficients
  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    std::cerr << "Polynomial interpolate node has no text near line " 
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate());
  }
  std::string sval = child->FIRST_TEXT_FN()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    std::cerr << "Polynomial interpolate node has invalid text near line "
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate()); 
  }

  std::vector<double> splits = split_string(sval);
  
  return std::shared_ptr<Interpolate>(new PolynomialInterpolate(splits));
}

std::shared_ptr<Interpolate> process_piecewise(const xmlpp::Element * child,
                                              int & ier)
{
  // The text is a intermoven list of (point, value) tuples
  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    std::cerr << "Polynomial interpolate node has no text near line " 
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate());
  }
  std::string sval = child->FIRST_TEXT_FN()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    std::cerr << "Polynomial interpolate node has invalid text near line "
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate()); 
  }

  std::vector<double> splits = split_string(sval);
  
  std::vector<double> points;
  std::vector<double> values;

  for (int i = 0; i < splits.size()/2; i++) {
    points.push_back(splits[i*2]);
    values.push_back(splits[i*2+1]);
  }

  return std::shared_ptr<Interpolate>(new PiecewiseLinearInterpolate(points, values));
}

std::shared_ptr<Interpolate> process_logpiecewise(const xmlpp::Element * child,
                                              int & ier)
{
  // The text is a intermoven list of (point, value) tuples
  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    std::cerr << "Log piecewise interpolate node has no text near line " 
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate());
  }
  std::string sval = child->FIRST_TEXT_FN()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    std::cerr << "Polynomial log piecewise interpolate node has invalid text near line "
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate()); 
  }

  std::vector<double> splits = split_string(sval);
  
  std::vector<double> points;
  std::vector<double> values;

  for (int i = 0; i < splits.size()/2; i++) {
    points.push_back(splits[i*2]);
    values.push_back(splits[i*2+1]);
  }

  return std::shared_ptr<Interpolate>(new PiecewiseLogLinearInterpolate(points, values));
}

std::shared_ptr<Interpolate> process_mts(const xmlpp::Element * child,
                                              int & ier)
{
  // The text is (V0, D, T0) -- probably fix later...
  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    std::cerr << "Polynomial interpolate node has no text near line " 
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate());
  }
  std::string sval = child->FIRST_TEXT_FN()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    std::cerr << "Polynomial interpolate node has invalid text near line "
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate()); 
  }

  std::vector<double> splits = split_string(sval);

  if (splits.size() != 3) {
    ier = BAD_TEXT;
    std::cerr << "MTS interpolate node does not contain 3 parameters near line "
        << child->get_line() << std::endl;
    return std::shared_ptr<Interpolate>(new InvalidInterpolate()); 
  }

  return std::shared_ptr<Interpolate>(new MTSShearInterpolate(
          splits[0], splits[1], splits[2]));
}

std::vector<double> split_string(std::string sval)
{
  std::vector<std::string> splits;
  std::stringstream ss(sval);
  std::string temp;
  while (ss >> temp) {
    splits.push_back(temp);
  }
  std::vector<double> value;
  for (auto it = splits.begin(); it != splits.end(); ++it) {
    value.push_back(std::stod(*it));
  }
  return value;
}

std::vector<std::shared_ptr<Interpolate>> 
vector_param(const xmlpp::Node * node, std::string name, int & ier)
{
  auto matches = node->get_children(name);
  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    std::cerr << "Multiple nodes of type " << name << " found near line " 
        << node->get_line() << std::endl;
    return {std::shared_ptr<Interpolate>(new InvalidInterpolate())};
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
    std::cerr << "Required node type " << name << " not found near line " 
        << node->get_line() << std::endl;
    return {std::shared_ptr<Interpolate>(new InvalidInterpolate())};
  }
  auto child = dynamic_cast<const xmlpp::Element*>(matches.front());
  
  // Two cases: text data -> vector of constants or multiple instances of
  // interpolate nodes
  auto inters = child->get_children("interpolate");
  if (inters.size() != 0) {
    std::vector<std::shared_ptr<Interpolate>> results;
    for (auto it = inters.begin(); it != inters.end(); ++it) {
      results.emplace_back(interpolate_node(*it, ier));
    }
    return results;
  }
  else if (child->has_child_text()) {
    std::string sval = child->FIRST_TEXT_FN()->get_content();
    if (sval == "") {
      ier = BAD_TEXT;
      std::cerr << "Vector parameter node near line " << child->get_line() << 
          " does not store text!" << std::endl;
      return {std::shared_ptr<Interpolate>(new InvalidInterpolate())};
    }
    std::vector<double> data = split_string(sval);
    std::vector<std::shared_ptr<Interpolate>> results;
    for (auto it = data.begin(); it != data.end(); ++it) {
      results.emplace_back(std::make_shared<ConstantInterpolate>(*it));
    }
    return results;
  }
  else {
    ier = BAD_TEXT;
    std::cerr << "Unable to process vector parameter near line " << 
        child->get_line() << std::endl;
    return {std::shared_ptr<Interpolate>(new InvalidInterpolate())};
  }
}

std::vector<double> 
vector_constant(const xmlpp::Node * node, std::string name, int & ier)
{
  auto matches = node->get_children(name);
  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    std::cerr << "Multiple nodes of type " << name << " found near line " 
        << node->get_line() << std::endl;
    return std::vector<double>();
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
    std::cerr << "Required node type " << name << " not found near line " 
        << node->get_line() << std::endl;
    return std::vector<double>();
  }
  auto child = dynamic_cast<const xmlpp::Element*>(matches.front());
  
  if (child->has_child_text()) {
    std::string sval = child->FIRST_TEXT_FN()->get_content();
    if (sval == "") {
      ier = BAD_TEXT;
      std::cerr << "Vector parameter node near line " << child->get_line() << 
          " does not store text!" << std::endl;
      return std::vector<double>();
    }
    return split_string(sval);
  }
  else {
    ier = BAD_TEXT;
    std::cerr << "Unable to process vector parameter near line " << 
        child->get_line() << std::endl;
    return std::vector<double>();
  }
}


} // namespace neml
