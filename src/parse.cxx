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

  // Find the named node
  std::stringstream ss;
  ss << "material[@name='" << mname << "']";
  auto fset = root->find(ss.str());

  // Default
  ier = SUCCESS;

  // Determine if we found anything
  if (fset.size() < 1) {
    ier = NODE_NOT_FOUND;
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  else if (fset.size() > 1) {
    ier = TOO_MANY_NODES;
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  else {
    return make_from_node(dynamic_cast<const xmlpp::Element*>(fset[0]), ier);
  }
}

std::shared_ptr<NEMLModel> make_from_node(const xmlpp::Element * node, int & ier)
{
  // Top level parse between types of materials
  std::string ft;
  if (not one_attribute(node, "type", ft, ier)) {
    return std::shared_ptr<NEMLModel>(nullptr);
  }

  if (ft == "smallstrain") {
    return process_smallstrain(node, ier);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<NEMLModel>(nullptr); 
  }
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
  // Top level parse between types of materials
  std::string ft;
  if (not one_attribute(node, "type", ft, ier)) {
    return std::shared_ptr<LinearElasticModel>(nullptr);
  }

  if (ft == "isotropic") {
    return process_isotropiclinearelastic(node, ier);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<LinearElasticModel>(nullptr); 
  } 
}

std::shared_ptr<LinearElasticModel> process_isotropiclinearelastic(
    const xmlpp::Element * node, int & ier)
{
  // Need a shear and a bulk modulus
  // Shear
  std::shared_ptr<ShearModulus> sm;
  xmlpp::Element * shear_node;
  if (not one_child(node, "shear", shear_node, ier)) {
    return std::shared_ptr<LinearElasticModel>(nullptr); 
  }
  sm = process_shearmodulus(shear_node, ier);
  if (ier != SUCCESS) return std::shared_ptr<LinearElasticModel>(nullptr);

  // Bulk
  std::shared_ptr<BulkModulus> bm;
  xmlpp::Element * bulk_node;
  if (not one_child(node, "bulk", bulk_node, ier)) {
    return std::shared_ptr<LinearElasticModel>(nullptr); 
  }
  bm = process_bulkmodulus(bulk_node, ier);
  if (ier != SUCCESS) return std::shared_ptr<LinearElasticModel>(nullptr);
  
  return std::shared_ptr<LinearElasticModel>(
      new IsotropicLinearElasticModel(sm,bm));
}

std::shared_ptr<ShearModulus> process_shearmodulus(
    const xmlpp::Element * node, int & ier)
{
  // Switch on the type of model
  std::string ft;
  if (not one_attribute(node, "type", ft, ier)) {
    return std::shared_ptr<ShearModulus>(nullptr);
  }

  if (ft == "constant") {
    return process_constantshearmodulus(node, ier);  
  }
  else if (ft == "polynomial") {
    return process_polynomialshearmodulus(node, ier);    
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<ShearModulus>(nullptr); 
  } 
}

std::shared_ptr<ShearModulus> process_constantshearmodulus(
    const xmlpp::Element * node, int & ier)
{
  // One scalar parameter in "modulus"
  double mod;
  if (not scalar_param(node, "modulus", mod, ier)) {
    return std::shared_ptr<ShearModulus>(nullptr);
  }

  return std::make_shared<ConstantShearModulus>(mod);
}

std::shared_ptr<ShearModulus> process_polynomialshearmodulus(
    const xmlpp::Element * node, int & ier)
{
  // One vector parameter in "coefs"
  std::vector<double> coefs;
  if (not vector_param(node, "coefs", coefs, ier)) {
    return std::shared_ptr<ShearModulus>(nullptr);
  }

  return std::make_shared<PolyShearModulus>(coefs);
}

std::shared_ptr<BulkModulus> process_bulkmodulus(
    const xmlpp::Element * node, int & ier)
{
  // Switch on the type of model
  std::string ft;
  if (not one_attribute(node, "type", ft, ier)) {
    return std::shared_ptr<BulkModulus>(nullptr);     
  }

  if (ft == "constant") {
    return process_constantbulkmodulus(node, ier);    

  }
  else if (ft == "polynomial") {
    return process_polynomialbulkmodulus(node, ier);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<BulkModulus>(nullptr); 
  } 
}

std::shared_ptr<BulkModulus> process_constantbulkmodulus(
    const xmlpp::Element * node, int & ier)
{
  // One scalar parameter in "modulus"
  double mod;
  if (not scalar_param(node, "modulus", mod, ier)) {
    return std::shared_ptr<BulkModulus>(nullptr);
  }

  return std::make_shared<ConstantBulkModulus>(mod);
}

std::shared_ptr<BulkModulus> process_polynomialbulkmodulus(
    const xmlpp::Element * node, int & ier)
{
  // One vector parameter in "coefs"
  std::vector<double> coefs;
  if (not vector_param(node, "coefs", coefs, ier)) {
    return std::shared_ptr<BulkModulus>(nullptr);
  }

  return std::make_shared<PolyBulkModulus>(coefs);
}


std::shared_ptr<RateIndependentFlowRule> process_independent(
    const xmlpp::Element * node, int & ier)
{
  // Should have a "rule" node of differing type
  xmlpp::Element * rule_node;
  if (not one_child(node, "rule", rule_node, ier)) {
    return std::shared_ptr<RateIndependentFlowRule>(nullptr);
  }
  return process_rirule(rule_node, ier);
}

std::shared_ptr<RateIndependentFlowRule> process_rirule(
    const xmlpp::Element * node, int & ier)
{
  // Types are "associative" or "nonassociative"
  std::string ft;
  if (not one_attribute(node, "type", ft, ier)) {
    return std::shared_ptr<RateIndependentFlowRule>(nullptr);
  }

  if (ft == "associative") {
    return process_riassociative(node, ier);
  }
  else if (ft == "nonassociative") {
    return process_rinonassociative(node, ier);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<RateIndependentFlowRule>(nullptr);
  }
}

std::shared_ptr<RateIndependentFlowRule> process_riassociative(
    const xmlpp::Element * node, int & ier)
{
  // Need a surface and a hardening model
  // Surface
  xmlpp::Element * surface_node;
  if (not one_child(node, "surface", surface_node, ier)) {
    return std::shared_ptr<RateIndependentFlowRule>(nullptr);
  }
  std::shared_ptr<YieldSurface> ys = process_surface(surface_node, ier);
  if (ier != SUCCESS) {
    return std::shared_ptr<RateIndependentFlowRule>(nullptr);
  }

  // Hardening model
  xmlpp::Element * hardening_node;
  if (not one_child(node, "hardening", hardening_node, ier)) {
    return std::shared_ptr<RateIndependentFlowRule>(nullptr);
  }
  std::shared_ptr<HardeningRule> hr = process_hardening(hardening_node, ier);
  if (ier != SUCCESS) {
    return std::shared_ptr<RateIndependentFlowRule>(nullptr);
  }

  return std::make_shared<RateIndependentAssociativeFlow>(ys, hr);
}

std::shared_ptr<YieldSurface> process_surface(
    const xmlpp::Element * node, int & ier)
{
  // Switch on type: isoj2, isokinj2
  std::string ft;
  if (not one_attribute(node, "type", ft, ier)) {
    return std::shared_ptr<YieldSurface>(nullptr);
  }
  
  if (ft == "isoj2") {
    return process_isoj2(node, ier);
  }
  else if (ft == "isokinj2") {
    return process_isokinj2(node, ier);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<YieldSurface>(nullptr);
  }
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
  std::string ft;
  if (not one_attribute(node, "type", ft, ier)) {
    return std::shared_ptr<HardeningRule>(nullptr);
  }
  
  if (ft == "isotropic") {
    return process_isotropic(node, ier);
  }
  else if (ft == "kinematic") {
    return process_kinematic(node, ier);
  }
  else if (ft == "combined") {
    return process_combined(node, ier);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<HardeningRule>(nullptr);
  }
}

std::shared_ptr<IsotropicHardeningRule> process_isotropic(
    const xmlpp::Element * node, int & ier)
{ 
  // Should have an "isotropic" tag with a type
  xmlpp::Element * isotropic_node;
  if (not one_child(node, "isotropic", isotropic_node, ier)) {
    return std::shared_ptr<IsotropicHardeningRule>(nullptr);
  }

  // Two options: linear or voce
  std::string ft;
  if (not one_attribute(isotropic_node, "type", ft, ier)) {
    return std::shared_ptr<IsotropicHardeningRule>(nullptr);
  }

  if (ft == "linear") {
    return process_linearisotropic(isotropic_node, ier);
  }
  else if (ft == "voce") {
    return process_voceisotropic(isotropic_node, ier);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<IsotropicHardeningRule>(nullptr);
  }
}

std::shared_ptr<IsotropicHardeningRule> process_linearisotropic(
    const xmlpp::Element * node, int & ier)
{
  // Two parameters: "yield" and "harden"
  double yield;
  if (not scalar_param(node, "yield", yield, ier)) {
    return std::shared_ptr<IsotropicHardeningRule>(nullptr);
  }

  double harden;
  if (not scalar_param(node, "harden", harden, ier)) {
    return std::shared_ptr<IsotropicHardeningRule>(nullptr);
  }

  return std::make_shared<LinearIsotropicHardeningRule>(yield, harden);
}

std::shared_ptr<IsotropicHardeningRule> process_voceisotropic(
    const xmlpp::Element * node, int & ier)
{
  // Three parameters: "yield", "r", and "d"
  double yield;
  if (not scalar_param(node, "yield", yield, ier)) {
    return std::shared_ptr<IsotropicHardeningRule>(nullptr);
  }

  double r;
  if (not scalar_param(node, "r", r, ier)) {
    return std::shared_ptr<IsotropicHardeningRule>(nullptr);
  }

  double d;
  if (not scalar_param(node, "d", d, ier)) {
    return std::shared_ptr<IsotropicHardeningRule>(nullptr);
  }

  return std::make_shared<VoceIsotropicHardeningRule>(yield, r, d);
}


std::shared_ptr<KinematicHardeningRule> process_kinematic(
    const xmlpp::Element * node, int & ier)
{
  // Should have a "kinematic" tag with a type
  xmlpp::Element * kinematic_node;
  if (not one_child(node, "kinematic", kinematic_node, ier)) {
    return std::shared_ptr<KinematicHardeningRule>(nullptr);
  }

  // Only one option: linear
  std::string ft;
  if (not one_attribute(kinematic_node, "type", ft, ier)) {
    return std::shared_ptr<KinematicHardeningRule>(nullptr);
  }
  
  if (ft == "linear") {
    return process_linearkinematic(kinematic_node, ier);
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<KinematicHardeningRule>(nullptr);
  }

}

std::shared_ptr<KinematicHardeningRule> process_linearkinematic(
    const xmlpp::Element * node, int & ier)
{
  // One parameter: "harden"
  double harden;
  if (not scalar_param(node, "harden", harden, ier)) {
    return std::shared_ptr<KinematicHardeningRule>(nullptr);
  }

  return std::make_shared<LinearKinematicHardeningRule>(harden);
}

std::shared_ptr<HardeningRule> process_combined(const xmlpp::Element * node,
                                                int & ier)
{
  // Recurse to pick up "isotropic" and "kinematic"

  // isotropic
  xmlpp::Element * isotropic_node;
  if (not one_child(node, "isotropic", isotropic_node, ier)) {
    return std::shared_ptr<HardeningRule>(nullptr);
  }
  std::shared_ptr<IsotropicHardeningRule> ir = process_isotropic(
      isotropic_node, ier);
  if (ier != SUCCESS) {
    return std::shared_ptr<HardeningRule>(nullptr);
  }

  // kinematic
  xmlpp::Element * kinematic_node;
  if (not one_child(node, "kinematic", kinematic_node, ier)) {
    return std::shared_ptr<HardeningRule>(nullptr);
  }
  std::shared_ptr<KinematicHardeningRule> kr = process_kinematic(
      kinematic_node, ier);
  if (ier != SUCCESS) {
    return std::shared_ptr<HardeningRule>(nullptr);
  }

  return std::make_shared<CombinedHardeningRule>(ir, kr);
}

std::shared_ptr<RateIndependentFlowRule> process_rinonassociative(
    const xmlpp::Element * node, int & ier)
{

}

std::shared_ptr<ViscoPlasticFlowRule> process_dependent(
    const xmlpp::Element * node, int & ier)
{
  // Should have a "rule" node of differing type
  xmlpp::Element * rule_node;
  if (not one_child(node, "rule", rule_node, ier)) {
    return std::shared_ptr<ViscoPlasticFlowRule>(nullptr);
  }
  return process_rdrule(rule_node, ier);
}

std::shared_ptr<ViscoPlasticFlowRule> process_rdrule(
    const xmlpp::Element * node, int & ier)
{
  // Types are "associative", "chaboche", or "yaguchigr91"
}


// Helpers
bool one_child(const xmlpp::Node * node, std::string name,
               xmlpp::Element * & child, int & ier)
{
  auto matches = node->get_children(name);

  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    return false;
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
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
    return false;
  }
  else {
    return true;
  }
}

bool scalar_param(const xmlpp::Node * node, std::string name,
                  double & value, int & ier)
{
  auto matches = node->get_children(name);
  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    return false;
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
    return false;
  }
  auto child = dynamic_cast<xmlpp::Element*>(matches.front());

  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    return false;
  }

  std::string sval = child->get_child_text()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    return false;
  }
  value = std::stod(sval);

  return true;
}

bool vector_param(const xmlpp::Node * node, std::string name,
                  std::vector<double> & value, int & ier)
{
  auto matches = node->get_children(name);
  if (matches.size() > 1) {
    ier = TOO_MANY_NODES;
    return false;
  }
  else if (matches.size() < 1) {
    ier = NODE_NOT_FOUND;
    return false;
  }
  auto child = dynamic_cast<xmlpp::Element*>(matches.front());

  if (not child->has_child_text()) {
    ier = BAD_TEXT;
    return false;
  }

  std::string sval = child->get_child_text()->get_content();
  if (sval == "") {
    ier = BAD_TEXT;
    return false;
  }

  std::vector<std::string> splits;
  std::stringstream ss(sval);
  std::string temp;
  while (ss >> temp) {
    splits.push_back(temp);
  }
  for (auto it = splits.begin(); it != splits.end(); ++it) {
    value.push_back(std::stod(*it));
  }

  return true;
}

} // namespace neml
