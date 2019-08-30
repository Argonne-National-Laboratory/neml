#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "objects.h"

#include <vector>
#include <memory>

namespace neml {

/// Base class for interpolation functions
//  This class defines a scalar interpolation function.
//  An implementation must also define the first derivative. 
class Interpolate: public NEMLObject {
 public:
  Interpolate();
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
class PolynomialInterpolate : public Interpolate {
 public:
  /// Input is the coefficients of the polynomial, from highest to lowest order
  PolynomialInterpolate(const std::vector<double> coefs);

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
class GenericPiecewiseInterpolate: public Interpolate {
 public:
  GenericPiecewiseInterpolate(std::vector<double> points,
                              std::vector<std::shared_ptr<Interpolate>> functions);

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
class PiecewiseLinearInterpolate: public Interpolate {
 public:
  /// Parameters are a list of x coordinates and a corresponding list of y 
  /// coordinates
  PiecewiseLinearInterpolate(const std::vector<double> points,
                             const std::vector<double> values);

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
class PiecewiseLogLinearInterpolate: public Interpolate {
 public:
  /// Similar to piecewise linear interpolation except the y coordinates are
  /// given as ln(y) and the interpolation is done in log space
  PiecewiseLogLinearInterpolate(const std::vector<double> points,
                             const std::vector<double> values);

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

/// A constant value
class ConstantInterpolate : public Interpolate {
 public:
  /// The parameter is the constant value!
  ConstantInterpolate(double v);

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
class ExpInterpolate : public Interpolate {
 public:
  /// The parameter is the constant value!
  ExpInterpolate(double A, double B);

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

/// The MTS shear modulus function proposed in the original paper
class MTSShearInterpolate : public Interpolate {
 public:
  /// Interpolation using the MTS model form f(x) = V0 - D / (exp(T0 / x) - 1)
  MTSShearInterpolate(double V0, double D, double T0);

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

/// A helper to make a vector of constant interpolates from a vector
std::vector<std::shared_ptr<Interpolate>> 
  make_vector(const std::vector<double> & iv);

/// A helper to evaluate a vector of interpolates
std::vector<double> eval_vector(
    const std::vector<std::shared_ptr<Interpolate>> & iv, double x);

/// A helper to evaluate the derivative of a vector of interpolates
std::vector<double> eval_deriv_vector(
    const std::vector<std::shared_ptr<Interpolate>> & iv, double x);

} // namespace neml


#endif // INTERPOLATE_H
