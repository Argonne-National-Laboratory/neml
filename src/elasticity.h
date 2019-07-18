#ifndef ELASTICITY_H
#define ELASTICITY_H

#include "objects.h"
#include "interpolate.h"
#include "math/tensors.h"
#include "math/rotations.h"

#include <memory>
#include <vector>
#include <string>
#include <set>

namespace neml {

/// Interface of all linear elastic models
//    Return properties as a function of temperature
class LinearElasticModel: public NEMLObject {
 public:
  LinearElasticModel();

  /// The stiffness tensor, in Mandel notation
  virtual int C(double T, double * const Cv) = 0;
  /// The compliance tensor, in Mandel notation
  virtual int S(double T, double * const Sv) = 0;

  /// The stiffness tensor in a tensor object
  SymSym C(double T);
  /// The compliance tensor in a tensor object
  SymSym S(double T);
  
  /// The rotated stiffness tensor in a tensor object
  SymSym C(double T, const Orientation & Q);
  /// The rotated compliance tensor in a tensor object
  SymSym S(double T, const Orientation & Q);

  /// An effective shear modulus
  virtual double G(double T);
  virtual double G(double T, const Orientation & Q, const Vector & b,
                   const Vector & n);

 private:
  void cache_orientation_(double T, const Orientation & Q);

 private:
  bool setup_;
  size_t hash_;
  SymSym C_;
  SymSym S_;

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
  virtual int C(double T, double * const Cv);
  /// Implement the compliance tensor
  virtual int S(double T, double * const Sv);

  /// The Young's modulus
  virtual double E(double T) const;
  /// Poisson's ratio
  virtual double nu(double T) const;
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

class CubicLinearElasticModel: public LinearElasticModel {
 public:
  CubicLinearElasticModel(std::shared_ptr<Interpolate> m1,
                          std::shared_ptr<Interpolate> m2,
                          std::shared_ptr<Interpolate> m3,
                          std::string method);

  /// The string type for the object system
  static std::string type();
  /// Setup default parameters for the object system
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();
  
  /// Implement the stiffness tensor
  virtual int C(double T, double * const Cv);
  /// Implement the compliance tensor
  virtual int S(double T, double * const Sv);

 private:
  void get_components_(double T, double & C1, double & C2, double & C3) const;


 private:
  std::shared_ptr<Interpolate> M1_, M2_, M3_;
  std::string method_;
};

static Register<CubicLinearElasticModel> regCubicLinearElasticModel;

} // namespace neml


#endif // ELASTICITY
