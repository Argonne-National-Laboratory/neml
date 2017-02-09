#include "interpolate.h"

#include "nemlmath.h"

#include <math.h>

namespace neml {

double Interpolate::operator()(double x) const
{
  return value(x);
}

PolynomialInterpolate::PolynomialInterpolate(const std::vector<double> coefs) :
    coefs_(coefs)
{

}

double PolynomialInterpolate::value(double x) const
{
  return polyval(&coefs_[0], coefs_.size(), x);
}

ConstantInterpolate::ConstantInterpolate(double v) :
    v_(v)
{

}

double ConstantInterpolate::value(double x) const
{
  return v_;
}

MTSShearInterpolate::MTSShearInterpolate(double V0, double D, double T0) :
    V0_(V0), D_(D), T0_(T0)
{
  
}

double MTSShearInterpolate::value(double x) const
{
  return V0_ - D_ / (exp(T0_ / x) - 1.0);
}

std::vector<std::shared_ptr<const Interpolate>> 
    make_constant_vector(const std::vector<double> & iv)
{
  std::vector<std::shared_ptr<const Interpolate>> vt;
  for (auto it = iv.begin(); it != iv.end(); ++it) {
    vt.emplace_back(std::make_shared<const ConstantInterpolate>(*it));
  }
  return vt;
}

std::vector<double> eval_vector(
    const std::vector<std::shared_ptr<const Interpolate>> & iv, double x)
{
  std::vector<double> vt;
  for (auto it = iv.begin(); it != iv.end(); ++it) {
    vt.push_back((*it)->value(x));
  }
  return vt;
}

} // namespace neml
