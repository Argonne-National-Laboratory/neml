#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>
#include <memory>

namespace neml {

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

class InvalidInterpolate : public Interpolate {
 public:
  InvalidInterpolate();
  virtual double value(double x) const;
  virtual double derivative(double x) const;
};

class PolynomialInterpolate : public Interpolate {
 public:
  PolynomialInterpolate(const std::vector<double> coefs);
  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> coefs_;
  std::vector<double> deriv_;

};

class PiecewiseLinearInterpolate: public Interpolate {
 public:
  PiecewiseLinearInterpolate(const std::vector<double> points,
                             const std::vector<double> values);
  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const std::vector<double> points_, values_;

};

class ConstantInterpolate : public Interpolate {
 public:
  ConstantInterpolate(double v);
  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double v_;

};

class MTSShearInterpolate : public Interpolate {
 public:
  MTSShearInterpolate(double V0, double D, double T0);
  virtual double value(double x) const;
  virtual double derivative(double x) const;

 private:
  const double V0_, D_, T0_;
};

std::vector<std::shared_ptr<Interpolate>> 
  make_vector(const std::vector<double> & iv);

std::vector<double> eval_vector(
    const std::vector<std::shared_ptr<Interpolate>> & iv, double x);

} // namespace neml


#endif // INTERPOLATE_H
