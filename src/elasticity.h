#ifndef ELASTICITY_H
#define ELASTICITY_H

#include <memory>
#include <vector>

namespace neml {

/// Shear modulus model: function of temperature
class ShearModulus {
 public:
  virtual double modulus(double T) const = 0;

};

/// Constant shear modulus
class ConstantShearModulus: public ShearModulus {
 public:
  ConstantShearModulus(double mu);
  virtual double modulus(double T) const;

  double mu() const;

 private:
  const double mu_;

};

/// Shear modulus varies as a general polynomial
class PolyShearModulus: public ShearModulus {
 public:
  PolyShearModulus(const std::vector<double> coefs);
  PolyShearModulus(int n, const double * const coefs);
  virtual double modulus(double T) const;

  // Getters
  int n() const;
  const std::vector<double> & coefs() const;

 private:
  const int n_;
  const std::vector<double> coefs_;
};

/// Bulk modulus model: function of temperature
class BulkModulus {
 public:
  virtual double modulus(double T) const = 0;
};

/// Constant bulk modulus
class ConstantBulkModulus {
 public:
  ConstantBulkModulus(double K);
  virtual double modulus(double T) const;

  double K() const;

 private:
  const double K_;
};

/// Bulk modulus varies as a general polynomial
class PolyBulkModulus: public BulkModulus {
 public:
  PolyBulkModulus(const std::vector<double> coefs);
  PolyBulkModulus(int n, const double * const coefs);
  virtual double modulus(double T) const;

  // Getters
  int n() const;
  const std::vector<double> & coefs() const;

 private:
  const int n_;
  const std::vector<double> coefs_;
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
