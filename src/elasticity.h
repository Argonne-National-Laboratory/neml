#ifndef ELASTICITY_H
#define ELASTICITY_H

#include "interpolate.h"

#include <memory>
#include <vector>

namespace neml {

/// Shear modulus model: function of temperature
class ShearModulus {
 public:
  ShearModulus::ShearModulus(double mu);
  ShearModulus::ShearModulus(std::shared_ptr<Interpolate> mu);
  virtual double modulus(double T) const;

 private:
  std::shared_ptr<Interpolate> mu_;

};

/// Bulk modulus model: function of temperature
class BulkModulus {
 public:
  BulkModulus::BulkModulus(double K);
  BulkModulus::BulkModulus(std::shared_ptr<Interpolate> K);
  virtual double modulus(double T) const;

 private:
  std::shared_ptr<Interpolate> K_;
};

/// Interface of all linear elastic models
//    Return properties as a function of temperature
class LinearElasticModel {
 public:
  virtual int C(double T, double * const Cv) const = 0;
  virtual int S(double T, double * const Sv) const = 0;
};

/// Isotropic shear modulus generating properties from shear and bulk models
class IsotropicLinearElasticModel: public LinearElasticModel {
  public:
   IsotropicLinearElasticModel(
       std::shared_ptr<ShearModulus> shear, 
       std::shared_ptr<BulkModulus> bulk);

   virtual int C(double T, double * const Cv) const;
   virtual int S(double T, double * const Sv) const;

  private:
   std::shared_ptr<ShearModulus> shear_;
   std::shared_ptr<BulkModulus> bulk_;
};

} // namespace neml


#endif // ELASTICITY
