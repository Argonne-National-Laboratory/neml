#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>
#include <memory>

namespace neml {

/// Base class for interpolation functions
//  This class defines a scalar interpolation function.
//  An implementation must also define the first derivative. 
class Interpolate {
 public:
  Interpolate();
  virtual double value(double x) const = 0;
  virtual double derivative(double x) const = 0;
  double operator()(double x) const;
  bool valid() const;

 protected:
  bool valid_;
};

/// A dummy interpolation.
class InvalidInterpolate : public Interpolate {
 public:
  InvalidInterpolate();
  virtual double value(double x) const;
  virtual double derivative(double x) const;
};

/// Simple polynomial interpolation
class PolynomialInterpolate : public Interpolate {
 public:
  PolynomialInterpolate(const std::vector<double> coefs);
  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> coefs_;
  std::vector<double> deriv_;

};

/// Piecewise linear interpolation
class PiecewiseLinearInterpolate: public Interpolate {
 public:
  PiecewiseLinearInterpolate(const std::vector<double> points,
                             const std::vector<double> values);
  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> points_, values_;

};

/// A constant value
class ConstantInterpolate : public Interpolate {
 public:
  ConstantInterpolate(double v);
  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double v_;

};

/// The MTS shear modulus function proposed in the original paper
class MTSShearInterpolate : public Interpolate {
 public:
  MTSShearInterpolate(double V0, double D, double T0);
  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double V0_, D_, T0_;
};

/// A helper to make a vector of constant interpolates from a vector
std::vector<std::shared_ptr<Interpolate>> 
  make_vector(const std::vector<double> & iv);

/// A helper to evaluate a vector of interpolates
std::vector<double> eval_vector(
    const std::vector<std::shared_ptr<Interpolate>> & iv, double x);

} // namespace neml


#endif // INTERPOLATE_H
