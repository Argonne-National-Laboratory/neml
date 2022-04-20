#ifndef ELASTICITY_H
#define ELASTICITY_H

#include "objects.h"
#include "interpolate.h"
#include "math/tensors.h"
#include "math/rotations.h"

#include "windows.h"

#include <memory>
#include <vector>
#include <string>
#include <set>

namespace neml {

/// Interface of all linear elastic models
//    Return properties as a function of temperature
class NEML_EXPORT LinearElasticModel: public NEMLObject {
 public:
  LinearElasticModel(ParameterSet & params);

  /// The stiffness tensor, in Mandel notation
  virtual void C(double T, double * const Cv) const = 0;
  /// The compliance tensor, in Mandel notation
  virtual void S(double T, double * const Sv) const = 0;

  /// The stiffness tensor in a tensor object
  SymSymR4 C(double T) const;
  /// The compliance tensor in a tensor object
  SymSymR4 S(double T) const;

  /// The rotated stiffness tensor in a tensor object
  SymSymR4 C(double T, const Orientation & Q) const;
  /// The rotated compliance tensor in a tensor object
  SymSymR4 S(double T, const Orientation & Q) const;

  /// An effective shear modulus
  virtual double G(double T) const;
  virtual double G(double T, const Orientation & Q, const Vector & b,
                   const Vector & n) const;

};

/// Isotropic shear modulus generating properties from shear and bulk models
class NEML_EXPORT IsotropicLinearElasticModel: public LinearElasticModel {
 public:
  /// See detailed documentation for how to initialize with elastic constants
  IsotropicLinearElasticModel(ParameterSet & params);

  /// The string type for the object system
  static std::string type();
  /// Setup default parameters for the object system
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();

  /// Implement the stiffness tensor
  virtual void C(double T, double * const Cv) const;
  /// Implement the compliance tensor
  virtual void S(double T, double * const Sv) const;

  /// The Young's modulus
  virtual double E(double T) const;
  /// Poisson's ratio
  virtual double nu(double T) const;
  /// The bulk modulus
  virtual double K(double T) const;

 private:
  void C_calc_(double G, double K, double * const Cv) const;
  void S_calc_(double G, double K, double * const Sv) const;

  void get_GK_(double T, double & G, double & K) const;

 private:
  std::shared_ptr<Interpolate> m1_, m2_;
  std::string m1_type_, m2_type_;
  const std::set<std::string> valid_types_ = {"bulk", "shear",
    "youngs", "poissons"};
};

static Register<IsotropicLinearElasticModel> regIsotropicLinearElasticModel;

class NEML_EXPORT CubicLinearElasticModel: public LinearElasticModel {
 public:
  CubicLinearElasticModel(ParameterSet & params);

  /// The string type for the object system
  static std::string type();
  /// Setup default parameters for the object system
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();

  /// Implement the stiffness tensor
  virtual void C(double T, double * const Cv) const;
  /// Implement the compliance tensor
  virtual void S(double T, double * const Sv) const;

 private:
  void get_components_(double T, double & C1, double & C2, double & C3) const;


 private:
  std::shared_ptr<Interpolate> M1_, M2_, M3_;
  std::string method_;
};

static Register<CubicLinearElasticModel> regCubicLinearElasticModel;

class NEML_EXPORT TransverseIsotropicLinearElasticModel: public LinearElasticModel {
 public:
  TransverseIsotropicLinearElasticModel(ParameterSet & params);

  /// The string type for the object system
  static std::string type();
  /// Setup default parameters for the object system
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();

  /// Implement the stiffness tensor
  virtual void C(double T, double * const Cv) const;
  /// Implement the compliance tensor
  virtual void S(double T, double * const Sv) const;

 private:
  void get_components_(double T, double & C11, double & C33, double & C12,
                       double & C13, double & C44) const;


 private:
  std::shared_ptr<Interpolate> m1_, m2_, m3_, m4_, m5_;
  std::string method_;
};

static Register<TransverseIsotropicLinearElasticModel> regTransverseIsotropicLinearElasticModel;

} // namespace neml


#endif // ELASTICITY
