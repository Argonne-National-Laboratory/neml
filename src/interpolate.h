#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>
#include <memory>

namespace neml {

class Interpolate {
 public:
  virtual double value(double x) const = 0;
  double operator()(double x) const;
};

class PolynomialInterpolate : public Interpolate {
 public:
  PolynomialInterpolate(const std::vector<double> coefs);
  virtual double value(double x) const;

 private:
  const std::vector<double> coefs_;

};

class ConstantInterpolate : public Interpolate {
 public:
  ConstantInterpolate(double v);
  virtual double value(double x) const;

 private:
  const double v_;

};

std::vector<std::shared_ptr<const Interpolate>> 
  make_constant_vector(const std::vector<double> & iv);

std::vector<double> eval_vector(
    const std::vector<std::shared_ptr<const Interpolate>> & iv, double x);

} // namespace neml


#endif // INTERPOLATE_H
