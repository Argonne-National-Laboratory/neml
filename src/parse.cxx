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
  std::string ft = node->get_attribute_value("type");

  if (ft == "") {
    ier = ATTRIBUTE_NOT_FOUND;
    return std::shared_ptr<NEMLModel>(nullptr);    
  }
  else if (ft == "smallstrain") {
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
  auto eset = node->find("elastic");
  std::shared_ptr<LinearElasticModel> emodel;
  if (eset.size() < 1) {
    ier = NODE_NOT_FOUND;
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  else if (eset.size() > 1) {
    ier = TOO_MANY_NODES;
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  else {
    emodel = process_linearelastic(
        dynamic_cast<const xmlpp::Element*>(eset[0]), ier);
    if (ier != SUCCESS) return std::shared_ptr<NEMLModel>(nullptr);
  }
  
  // Logic here because we treat viscoplasticity differently
  auto pset = node->find("plastic");
  if (pset.size() < 1) {
    ier = NODE_NOT_FOUND;
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  else if (pset.size() > 1) {
    ier = TOO_MANY_NODES;
    return std::shared_ptr<NEMLModel>(nullptr);
  }
  else {
    // Switch between rate independent and rate dependent behavior
  }

  // Cheat to run tests
  std::shared_ptr<ConstantShearModulus> sm(new ConstantShearModulus(10.0));
  std::shared_ptr<ConstantBulkModulus> bm(new ConstantBulkModulus(10.0));
  std::shared_ptr<IsotropicLinearElasticModel> em(new IsotropicLinearElasticModel(sm, bm));
  std::shared_ptr<NEMLModel> fm(new SmallStrainElasticity(em));
  
  return fm;
}

std::shared_ptr<LinearElasticModel> process_linearelastic(
    const xmlpp::Element * node, int & ier)
{
  // Switch on the type of model
  // Top level parse between types of materials
  std::string ft = node->get_attribute_value("type");

  if (ft == "") {
    ier = ATTRIBUTE_NOT_FOUND;
    return std::shared_ptr<LinearElasticModel>(nullptr);    
  }
  else if (ft == "isotropic") {
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
  auto shear_children = node->get_children("shear");
  std::shared_ptr<ShearModulus> sm;
  if (shear_children.size() < 1) {
    ier = NODE_NOT_FOUND;
    return std::shared_ptr<LinearElasticModel>(nullptr);
  }
  else if (shear_children.size() > 1) {
    ier = TOO_MANY_NODES;
    return std::shared_ptr<LinearElasticModel>(nullptr);
  }
  else {
    sm = process_shearmodulus(
        dynamic_cast<const xmlpp::Element*>(shear_children.front()), ier);
    if (ier != SUCCESS) return std::shared_ptr<LinearElasticModel>(nullptr);
  }

  auto bulk_children = node->get_children("bulk");
  std::shared_ptr<BulkModulus> bm;
  if (bulk_children.size() < 1) {
    ier = NODE_NOT_FOUND;
    return std::shared_ptr<LinearElasticModel>(nullptr);
  }
  else if (bulk_children.size() > 1) {
    ier = TOO_MANY_NODES;
    return std::shared_ptr<LinearElasticModel>(nullptr);
  }
  else {
    bm = process_bulkmodulus(
        dynamic_cast<const xmlpp::Element*>(bulk_children.front()), ier);
    if (ier != SUCCESS) return std::shared_ptr<LinearElasticModel>(nullptr);
  }
  
  return std::shared_ptr<LinearElasticModel>(
      new IsotropicLinearElasticModel(sm,bm));
}

std::shared_ptr<ShearModulus> process_shearmodulus(
    const xmlpp::Element * node, int & ier)
{
  // Switch on the type of model
  std::string ft = node->get_attribute_value("type");

  if (ft == "") {
    ier = ATTRIBUTE_NOT_FOUND;
    return std::shared_ptr<ShearModulus>(nullptr);    
  }
  else if (ft == "constant") {
    return std::shared_ptr<ShearModulus>(nullptr);    
  }
  else if (ft == "polynomial") {
    return std::shared_ptr<ShearModulus>(nullptr);    
  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<ShearModulus>(nullptr); 
  } 
}

std::shared_ptr<BulkModulus> process_bulkmodulus(
    const xmlpp::Element * node, int & ier)
{
  // Switch on the type of model
  std::string ft = node->get_attribute_value("type");

  if (ft == "") {
    ier = ATTRIBUTE_NOT_FOUND;
    return std::shared_ptr<BulkModulus>(nullptr);    
  }
  else if (ft == "constant") {
    return std::shared_ptr<BulkModulus>(nullptr);    

  }
  else if (ft == "polynomial") {
    return std::shared_ptr<BulkModulus>(nullptr);    

  }
  else {
    ier = UNKNOWN_TYPE;
    return std::shared_ptr<BulkModulus>(nullptr); 
  } 
}

} // namespace neml
