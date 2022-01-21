#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "objects.h"

#include "windows.h"

#include <vector>
#include <memory>

#include <iostream>

namespace neml {

/// Base class for interpolation functions
//  This class defines a scalar interpolation function.
//  An implementation must also define the first derivative.
class NEML_EXPORT Interpolate: public NEMLObject {
 public:
  Interpolate(ParameterSet & params);
  /// Returns the value of the function
  virtual double value(double x) const = 0;
  /// Returns the derivative of the function
  virtual double derivative(double x) const = 0;
  /// Nice wrapper for function call syntax
  double operator()(double x) const;
  /// Is the interpolate valid?
  bool valid() const;

 protected:
  bool valid_;
};

/// Simple polynomial interpolation
class NEML_EXPORT PolynomialInterpolate : public Interpolate {
 public:
  /// Input is the coefficients of the polynomial, from highest to lowest order
  PolynomialInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> coefs_;
  std::vector<double> deriv_;
};

static Register<PolynomialInterpolate> regPolynomialInterpolate;

/// Generic piecewise interpolation
class NEML_EXPORT GenericPiecewiseInterpolate: public Interpolate {
 public:
  GenericPiecewiseInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> points_;
  const std::vector<std::shared_ptr<Interpolate>> functions_;
};

static Register<GenericPiecewiseInterpolate> regGenericPiecewiseInterpolate;

/// Piecewise linear interpolation
class NEML_EXPORT PiecewiseLinearInterpolate: public Interpolate {
 public:
  /// Parameters are a list of x coordinates and a corresponding list of y
  /// coordinates
  PiecewiseLinearInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> points_, values_;
};

static Register<PiecewiseLinearInterpolate> regPiecewiseLinearInterpolate;

/// Piecewise loglinear interpolation
class NEML_EXPORT PiecewiseLogLinearInterpolate: public Interpolate {
 public:
  /// Similar to piecewise linear interpolation except the y coordinates are
  /// given as ln(y) and the interpolation is done in log space
  PiecewiseLogLinearInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> points_;
  std::vector<double> values_;
};

static Register<PiecewiseLogLinearInterpolate> regPiecewiseLogLinearInterpolate;

/// Piecewise semiloglinear interpolation
class NEML_EXPORT PiecewiseSemiLogXLinearInterpolate: public Interpolate {
 public:
  /// Similar to piecewise linear interpolation except the interpolation is done
  /// in log space
  PiecewiseSemiLogXLinearInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> points_;
  std::vector<double> values_;
};

static Register<PiecewiseSemiLogXLinearInterpolate> regPiecewiseSemiLogXLinearInterpolate;

/// A constant value
class NEML_EXPORT ConstantInterpolate : public Interpolate {
 public:
  /// The parameter is the constant value!
  ConstantInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double v_;
};

static Register<ConstantInterpolate> regConstantInterpolate;

/// A*exp(B/x)
class NEML_EXPORT ExpInterpolate : public Interpolate {
 public:
  /// A and B parameters for formula above
  ExpInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double A_, B_;
};

static Register<ExpInterpolate> regExpInterpolate;

/// A*x**B
class NEML_EXPORT PowerLawInterpolate : public Interpolate {
 public:
  /// A and B for the exponential formula above
  PowerLawInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double A_, B_;
};

static Register<PowerLawInterpolate> regPowerLawInterpolate;

/// The MTS shear modulus function proposed in the original paper
class NEML_EXPORT MTSShearInterpolate : public Interpolate {
 public:
  /// Interpolation using the MTS model form f(x) = V0 - D / (exp(T0 / x) - 1)
  MTSShearInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double V0_, D_, T0_;
};

static Register<MTSShearInterpolate> regMTSShearInterpolate;


/// The Mechanical Threshold Stress scaling
class NEML_EXPORT MTSInterpolate : public Interpolate {
 public:
  /// Interpolation using the MTS model
  MTSInterpolate(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Create parameters for the object system
  static ParameterSet parameters();
  /// Create object from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double tau0_, g0_, q_, p_, k_, b_;
  const std::shared_ptr<Interpolate> mu_;
};

static Register<MTSInterpolate> regMTSInterpolate;

/// A helper to make a vector of constant interpolates from a vector
NEML_EXPORT std::vector<std::shared_ptr<Interpolate>>
  make_vector(const std::vector<double> & iv);

/// A helper to evaluate a vector of interpolates
NEML_EXPORT std::vector<double> eval_vector(
    const std::vector<std::shared_ptr<Interpolate>> & iv, double x);

/// A helper to evaluate the derivative of a vector of interpolates
NEML_EXPORT std::vector<double> eval_deriv_vector(
    const std::vector<std::shared_ptr<Interpolate>> & iv, double x);

/// A helper to make a constant interpolate from a double
NEML_EXPORT std::shared_ptr<ConstantInterpolate> make_constant(double v);

/// A helper to make a constant interpolate from a double
NEML_EXPORT std::unique_ptr<ConstantInterpolate> make_constant_unique(double v);

/// A helper to make a piecewiselinear interpolate
NEML_EXPORT std::shared_ptr<PiecewiseLinearInterpolate> make_piecewise(
    std::vector<double> points, std::vector<double> values);

} // namespace neml


#endif // INTERPOLATE_H
