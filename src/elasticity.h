#ifndef ELASTICITY_H
#define ELASTICITY_H

#include "objects.h"
#include "interpolate.h"

#include <memory>
#include <vector>
#include <string>
#include <set>

namespace neml {

/// Interface of all linear elastic models
//    Return properties as a function of temperature
class LinearElasticModel: public NEMLObject {
 public:
  /// The stiffness tensor, in Mandel notation
  virtual int C(double T, double * const Cv) const = 0;
  /// The compliance tensor, in Mandel notation
  virtual int S(double T, double * const Sv) const = 0;
  
  /// The Young's modulus
  virtual double E(double T) const = 0;
  /// Poisson's ratio
  virtual double nu(double T) const = 0;
  /// The shear modulus
  virtual double G(double T) const = 0;
  /// The bulk modulus
  virtual double K(double T) const = 0;
};

/// Isotropic shear modulus generating properties from shear and bulk models
class IsotropicLinearElasticModel: public LinearElasticModel {
 public:
  /// See detailed documentation for how to initialize with elastic constants
  IsotropicLinearElasticModel(
      std::shared_ptr<Interpolate> m1, 
      std::string m1_type,
      std::shared_ptr<Interpolate> m2,
      std::string m2_type);
  
  /// The string type for the object system
  static std::string type();
  /// Setup default parameters for the object system
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();
  
  /// Implement the stiffness tensor
  virtual int C(double T, double * const Cv) const;
  /// Implement the compliance tensor
  virtual int S(double T, double * const Sv) const;
  
  /// The Young's modulus
  virtual double E(double T) const;
  /// Poisson's ratio
  virtual double nu(double T) const;
  /// The shear modulus
  virtual double G(double T) const;
  /// The bulk modulus
  virtual double K(double T) const;
  
 private:
  int C_calc_(double G, double K, double * const Cv) const;
  int S_calc_(double G, double K, double * const Sv) const;

  void get_GK_(double T, double & G, double & K) const;
  
 private:
  std::shared_ptr<Interpolate> m1_, m2_;
  std::string m1_type_, m2_type_;
  const std::set<std::string> valid_types_ = {"bulk", "shear", 
    "youngs", "poissons"};
};

static Register<IsotropicLinearElasticModel> regIsotropicLinearElasticModel;

} // namespace neml


#endif // ELASTICITY
