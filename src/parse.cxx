#include "parse.h"

#include "nemlerror.h"

#include <sstream>
#include <iostream>

namespace neml {

std::shared_ptr<NEMLModel> parse_xml(std::string fname, std::string mname,
                                     int & ier)
{
  // Parse the XML file
  xmlpp::DomParser parser;
  parser.parse_file(fname);

  // Grab the root node
  const auto root = parser.get_document()->get_root_node();

  // Default
  ier = SUCCESS;

  // Dispatch to the right node
  return find_and_dispatch(root, "material", "name", mname, 
                           &make_from_node, ier);
}

std::shared_ptr<NEMLModel> make_from_node(const xmlpp::Element * node, int & ier)
{
  return dispatch_attribute<NEMLModel>(node, "type",
                                          {"smallstrain"},
                                          {&process_smallstrain},
                                          ier);
}

std::shared_ptr<NEMLModel> process_smallstrain(const xmlpp::Element * node, int & ier)
{
  // Need a elastic node and a plastic node
  // Elastic
  std::shared_ptr<LinearElasticModel> emodel;
  xmlpp::Element * elastic_node;
  if (not one_child(node, "elastic", elastic_node, ier)) {
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  emodel = process_linearelastic(elastic_node, ier);
  if (ier != SUCCESS) return std::shared_ptr<NEMLModel>(nullptr);
  
  // Logic here because we treat viscoplasticity differently
  xmlpp::Element * plastic_node;
  if (not one_child(node, "plastic", plastic_node, ier)) {
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  
  // Select then on independent/dependent
  std::string ft;
  if (not one_attribute(plastic_node, "type", ft, ier)) {
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  if (ft == "independent") {
    std::shared_ptr<RateIndependentFlowRule> fr =  process_independent(
        plastic_node, ier);
    
    return std::make_shared<SmallStrainRateIndependentPlasticity>(emodel, fr);
  }
  else if (ft == "dependent") {
    std::shared_ptr<ViscoPlasticFlowRule> fr = process_dependent(plastic_node, ier);
    std::shared_ptr<TVPFlowRule> gfr = std::make_shared<TVPFlowRule>(emodel, fr);
    return std::make_shared<GeneralIntegrator>(gfr);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<NEMLModel>(nullptr);
  }
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
  // Need a shear and a bulk modulus
  // Shear
  std::shared_ptr<ShearModulus> sm = dispatch_node(node, "shear", 
                                                   &process_shearmodulus, ier);

  // Bulk
  std::shared_ptr<BulkModulus> bm = dispatch_node(node, "bulk", 
                                                   &process_bulkmodulus, ier);
  
  return std::shared_ptr<LinearElasticModel>(
      new IsotropicLinearElasticModel(sm,bm));
}

std::shared_ptr<ShearModulus> process_shearmodulus(
    const xmlpp::Element * node, int & ier)
{
  // Switch on the type of model
  return dispatch_attribute<ShearModulus>(node, "type",
                                          {"constant", "polynomial"},
                                          {&process_constantshearmodulus, &process_polynomialshearmodulus},
                                          ier);
}

std::shared_ptr<ShearModulus> process_constantshearmodulus(
    const xmlpp::Element * node, int & ier)
{
  // One scalar parameter in "modulus"
  double mod = scalar_param(node, "modulus", ier);

  return std::make_shared<ConstantShearModulus>(mod);
}

std::shared_ptr<ShearModulus> process_polynomialshearmodulus(
    const xmlpp::Element * node, int & ier)
{
  // One vector parameter in "coefs"
  std::vector<double> coefs = vector_param(node, "coefs", ier);

  return std::make_shared<PolyShearModulus>(coefs);
}

std::shared_ptr<BulkModulus> process_bulkmodulus(
    const xmlpp::Element * node, int & ier)
{
  // Switch on the type of model
  return dispatch_attribute<BulkModulus>(node, "type",
                                          {"constant", "polynomial"},
                                          {&process_constantbulkmodulus, &process_polynomialbulkmodulus},
                                          ier);
}

std::shared_ptr<BulkModulus> process_constantbulkmodulus(
    const xmlpp::Element * node, int & ier)
{
  // One scalar parameter in "modulus"
  double mod = scalar_param(node, "modulus", ier);

  return std::make_shared<ConstantBulkModulus>(mod);
}

std::shared_ptr<BulkModulus> process_polynomialbulkmodulus(
    const xmlpp::Element * node, int & ier)
{
  // One vector parameter in "coefs"
  std::vector<double> coefs = vector_param(node, "coefs", ier);

  return std::make_shared<PolyBulkModulus>(coefs);
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
                                          {"isoj2", "isokinj2"},
                                          {&process_isoj2, &process_isokinj2},
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
                                          {"linear", "voce"},
                                          {&process_linearisotropic, &process_voceisotropic},
                                          ier);
}

std::shared_ptr<HardeningRule> process_linearisotropic(
    const xmlpp::Element * node, int & ier)
{
  // Two parameters: "yield" and "harden"
  double yield = scalar_param(node, "yield", ier);

  double harden = scalar_param(node, "harden", ier);

  return std::make_shared<LinearIsotropicHardeningRule>(yield, harden);
}

std::shared_ptr<HardeningRule> process_voceisotropic(
    const xmlpp::Element * node, int & ier)
{
  // Three parameters: "yield", "r", and "d"
  double yield = scalar_param(node, "yield", ier);

  double r = scalar_param(node, "r", ier);

  double d = scalar_param(node, "d", ier);

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
  double harden = scalar_param(node, "harden", ier);

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
  std::vector<double> c = vector_param(node, "c", ier);
  
  // vector of gamma models
  std::vector< std::shared_ptr<GammaModel> > gvec = 
      dispatch_vector_models<GammaModel>(
          node, "gammamodels", "gammamodel", &process_gammamodel, ier);
  
  // Create our object
  return std::make_shared<Chaboche>(cir, c, gvec);
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
  double g = scalar_param(node, "g", ier);

  return std::make_shared<ConstantGamma>(g);
}

std::shared_ptr<GammaModel> process_saturating_gamma(
    const xmlpp::Element * node, int & ier)
{
  // Three scalar parameters: gs, g0, beta
  double gs = scalar_param(node, "gs", ier);
  double g0 = scalar_param(node, "g0", ier);
  double beta = scalar_param(node, "beta", ier);

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
  double n = scalar_param(node, "n", ier);

  return std::make_shared<ChabocheFlowRule>(ys, nsr, fm, n);

}

std::shared_ptr<FluidityModel> process_fluidity(
    const xmlpp::Element * node, int & ier)
{
  // Right now only a constant fluidity option
  return dispatch_attribute<FluidityModel>(
      node, "type", {"constant"}, {&process_constant_fluidity}, ier);
}

std::shared_ptr<FluidityModel> process_constant_fluidity(
    const xmlpp::Element * node, int & ier)
{
  // One property, eta
  double eta = scalar_param(node, "eta", ier);

  return std::make_shared<ConstantFluidity>(eta);
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
  double eta = scalar_param(node, "eta", ier);

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
  double n = scalar_param(node, "n", ier);

  return std::make_shared<GPowerLaw>(n);
}



// Helpers
bool one_child(const xmlpp::Node * node, std::string name,
               xmlpp::Element * & child, int & ier)
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
    child = dynamic_cast<xmlpp::Element*>(matches.front());
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

double scalar_param(const xmlpp::Node * node, std::string name, int & ier)
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
  auto child = dynamic_cast<xmlpp::Element*>(matches.front());

  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    std::cerr << "Node " << name << " near line " 
        << node->get_line() << " does not appear to store data" << std::endl;
    return 0.0;
  }

  std::string sval = child->get_child_text()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    std::cerr << "Node " << name << " near line " 
        << node->get_line() << " does not appear to store data" << std::endl;
    return 0.0;
  }
  return std::stod(sval);
}

std::vector<double> vector_param(const xmlpp::Node * node, std::string name,
                            int & ier)
{
  auto matches = node->get_children(name);
  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    std::cerr << "Multiple nodes of type " << name << " found near line " 
        << node->get_line() << std::endl;
    return {0.0};
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
    std::cerr << "Required node type " << name << " not found near line " 
        << node->get_line() << std::endl;
    return {0.0};
  }
  auto child = dynamic_cast<xmlpp::Element*>(matches.front());

  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    std::cerr << "Node " << name << " near line " 
        << node->get_line() << " does not appear to store data" << std::endl;
    return {0.0};
  }

  std::string sval = child->get_child_text()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    std::cerr << "Node " << name << " near line " 
        << node->get_line() << " does not appear to store data" << std::endl;
    return {0.0};
  }

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

} // namespace neml
