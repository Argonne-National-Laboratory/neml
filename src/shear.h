#ifndef SHEAR_H
#define SHEAR_H

namespace neml {

/// Base class for shear modulus models
//    These give the shear modulus as a function of temperature
class ShearModulus {
 public:
  ShearModulus();
  virtual ~ShearModulus();

  virtual double modulus(double T_np1) = 0;

};

/// Constant, temperature independent shear modulus
class ConstantShearModulus: public ShearModulus {
 public:
  ConstantShearModulus(double mu);
  virtual ~ConstantShearModulus();

  virtual double modulus(double T_np1);

  // Getters
  double mu() const;

 private:
  const double mu_;

};


} // namespace neml

#endif // SHEAR_H
