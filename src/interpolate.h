#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>

namespace neml {

class Interpolate {
 public:
  virtual double value(double x) = 0;
  double operator()(double x);
};

class PolynomialInterpolate : public Interpolate {
 public:
  PolynomialInterpolate(const std::vector<double> coefs);
  virtual double value(double x);

 private:
  const std::vector<double> coefs_;

};

class ConstantInterpolate : public Interpolate {
 public:
  ConstantInterpolate(double v);
  virtual double value(double x);

 private:
  const double v_;

};

} // namespace neml


#endif // INTERPOLATE_H
