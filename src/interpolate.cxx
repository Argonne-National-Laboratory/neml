#include "interpolate.h"

#include "nemlmath.h"

namespace neml {

double Interpolate::operator()(double x)
{
  return value(x);
}

PolynomialInterpolate::PolynomialInterpolate(const std::vector<double> coefs) :
    coefs_(coefs)
{

}

double PolynomialInterpolate::value(double x)
{
  return polyval(&coefs_[0], coefs_.size(), x);
}

ConstantInterpolate::ConstantInterpolate(double v) :
    v_(v)
{

}

double ConstantInterpolate::value(double x)
{
  return v_;
}


} // namespace neml
