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
  virtual int C(double T, double * const Cv) const = 0;
  virtual int S(double T, double * const Sv) const = 0;

  virtual double E(double T) const = 0;
  virtual double nu(double T) const = 0;
  virtual double G(double T) const = 0;
  virtual double K(double T) const = 0;

  virtual bool valid() const;
};

/// Isotropic shear modulus generating properties from shear and bulk models
class IsotropicLinearElasticModel: public LinearElasticModel {
 public:
  IsotropicLinearElasticModel(
      std::shared_ptr<Interpolate> m1, 
      std::string m1_type,
      std::shared_ptr<Interpolate> m2,
      std::string m2_type);

  static std::string type();
  static std::shared_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual int C(double T, double * const Cv) const;
  virtual int S(double T, double * const Sv) const;

  virtual double E(double T) const;
  virtual double nu(double T) const;
  virtual double G(double T) const;
  virtual double K(double T) const;

  virtual bool valid() const;

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

class BlankElasticModel: public LinearElasticModel {
 public:
  BlankElasticModel();

  static std::string type();
  static std::shared_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual int C(double T, double * const Cv) const;
  virtual int S(double T, double * const Sv) const;

  virtual double E(double T) const;
  virtual double nu(double T) const;
  virtual double G(double T) const;
  virtual double K(double T) const;

};

static Register<BlankElasticModel> regBlankElasticModel;

} // namespace neml


#endif // ELASTICITY
