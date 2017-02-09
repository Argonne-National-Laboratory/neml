#ifndef ELASTICITY_H
#define ELASTICITY_H

#include "interpolate.h"

#include <memory>
#include <vector>

namespace neml {

/// Modulus parent class
class Modulus {
 public:
  Modulus(double x);
  Modulus(std::shared_ptr<Interpolate> x);
  virtual double modulus(double T) const;

 protected:
  std::shared_ptr<Interpolate> x_;
};

/// Shear modulus model: function of temperature
class ShearModulus: public Modulus {
 public:
  ShearModulus(double mu);
  ShearModulus(std::shared_ptr<Interpolate> mu);
};

/// Bulk modulus model: function of temperature
class BulkModulus: public Modulus {
 public:
  BulkModulus(double K);
  BulkModulus(std::shared_ptr<Interpolate> K);
};

/// Youngs modulus model: function of temperature
class YoungsModulus: public Modulus {
 public:
  YoungsModulus(double E);
  YoungsModulus(std::shared_ptr<Interpolate> E);
};

/// Poissons ratio model: function of temperature
class PoissonsRatio: public Modulus {
 public:
  PoissonsRatio(double nu);
  PoissonsRatio(std::shared_ptr<Interpolate> nu);
};


/// Interface of all linear elastic models
//    Return properties as a function of temperature
class LinearElasticModel {
 public:
  virtual int C(double T, double * const Cv) const = 0;
  virtual int S(double T, double * const Sv) const = 0;

   virtual double E(double T) const = 0;
   virtual double nu(double T) const = 0;
   virtual double G(double T) const = 0;
   virtual double K(double T) const = 0;
};

/// Isotropic shear modulus generating properties from shear and bulk models
class IsotropicLinearElasticModel: public LinearElasticModel {
  public:
   IsotropicLinearElasticModel(
       std::shared_ptr<ShearModulus> shear, 
       std::shared_ptr<BulkModulus> bulk);
   IsotropicLinearElasticModel(
       std::shared_ptr<YoungsModulus> youngs,
       std::shared_ptr<PoissonsRatio> poissons);

   virtual int C(double T, double * const Cv) const;
   virtual int S(double T, double * const Sv) const;

   virtual double E(double T) const;
   virtual double nu(double T) const;
   virtual double G(double T) const;
   virtual double K(double T) const;

  private:
   int C_calc_(double G, double K, double * const Cv) const;
   int S_calc_(double G, double K, double * const Sv) const;

   void get_GK_(double T, double & G, double & K) const;
  
   enum modulii { GK, Ev };
  
  private:
   std::shared_ptr<Modulus> a_, b_;
   const modulii have_modulii_;
};

} // namespace neml


#endif // ELASTICITY
