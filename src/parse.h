#ifndef PARSE_H
#define PARSE_H

#include "neml.h"

#include <string>
#include <vector>

#include <libxml++/libxml++.h>

namespace neml {

/// Main entry to parse model from xml file
std::shared_ptr<NEMLModel> parse_xml(std::string fname, std::string mname, 
                                     int & ier);

/// Setup a model from a root node
std::shared_ptr<NEMLModel> make_from_node(const xmlpp::Element * node, int & ier);

/// Setup a small strain model
std::shared_ptr<NEMLModel> process_smallstrain(const xmlpp::Element * node, int & ier);

/// Setup a linear elasticity model
std::shared_ptr<LinearElasticModel> process_linearelastic(const xmlpp::Element * node, int & ier);

/// Setup an isotropic linear elastic model
std::shared_ptr<LinearElasticModel> process_isotropiclinearelastic(const xmlpp::Element * node, int & ier);

/// Setup shear modulus
std::shared_ptr<ShearModulus> process_shearmodulus(const xmlpp::Element * node, int & ier);

/// Constant shear
std::shared_ptr<ShearModulus> process_constantshearmodulus(const xmlpp::Element * node, int & ier);

/// Polynomial shear
std::shared_ptr<ShearModulus> process_polynomialshearmodulus(const xmlpp::Element * node, int & ier);

/// Setup bulk modulus
std::shared_ptr<BulkModulus> process_bulkmodulus(const xmlpp::Element * node, int & ier);

/// Constant bulk
std::shared_ptr<BulkModulus> process_constantbulkmodulus(const xmlpp::Element * node, int & ier);

/// Polynomial bulk
std::shared_ptr<BulkModulus> process_polynomialbulkmodulus(const xmlpp::Element * node, int & ier);

/// Rate independent plasticity processing
std::shared_ptr<RateIndependentFlowRule> process_independent(const xmlpp::Element * node, int & ier);

/// Specific RI model
std::shared_ptr<RateIndependentFlowRule> process_rirule(const xmlpp::Element * node, int & ier);

/// Associative RI
std::shared_ptr<RateIndependentFlowRule> process_riassociative(const xmlpp::Element * node, int & ier);

/// Yield surfaces
std::shared_ptr<YieldSurface> process_surface(const xmlpp::Element * node, int & ier);

/// Isotropic J2
std::shared_ptr<YieldSurface> process_isoj2(const xmlpp::Element * node, int & ier);

/// Kinematic J2
std::shared_ptr<YieldSurface> process_isokinj2(const xmlpp::Element * node, int & ier);

/// Hardening rules
std::shared_ptr<HardeningRule> process_hardening(const xmlpp::Element * node, int & ier);

/// Isotropic hardening
std::shared_ptr<IsotropicHardeningRule> process_isotropic(const xmlpp::Element * node, int & ier);

/// Linear isotropic hardening
std::shared_ptr<IsotropicHardeningRule> process_linearisotropic(const xmlpp::Element * node, int & ier);

/// Voce isotropic hardening
std::shared_ptr<IsotropicHardeningRule> process_voceisotropic(const xmlpp::Element * node, int & ier);

/// Kinematic hardening
std::shared_ptr<KinematicHardeningRule> process_kinematic(const xmlpp::Element * node, int & ier);

/// Linear kinematic hardening
std::shared_ptr<KinematicHardeningRule> process_linearkinematic(const xmlpp::Element * node, int & ier);

/// Combined hardening
std::shared_ptr<HardeningRule> process_combined(const xmlpp::Element * node, int & ier);

/// Nonassociative RI
std::shared_ptr<RateIndependentFlowRule> process_rinonassociative(const xmlpp::Element * node, int & ier);

/// Rate dependent plasticity processing
std::shared_ptr<ViscoPlasticFlowRule> process_dependent(const xmlpp::Element * node, int & ier);

/// Specific RD model
std::shared_ptr<ViscoPlasticFlowRule> process_rdrule(const xmlpp::Element * node, int & ier);

// Begin helpers

/// Return a single child node with given name, if not return relevant error
bool one_child(const xmlpp::Node * node, std::string name,
               xmlpp::Element * & child, int & ier);

/// Return the string value of an attribute with a given name, if not error
bool one_attribute(const xmlpp::Element * node, std::string name,
                   std::string & value, int & ier);

/// Return the scalar value of a named child node
bool scalar_param(const xmlpp::Node * node, std::string name,
                  double & value, int & ier);

/// Return the vector value of a named child node
bool vector_param(const xmlpp::Node * node, std::string name,
                  std::vector<double> & value, int & ier);

} // namespace neml

#endif // PARSE_H
