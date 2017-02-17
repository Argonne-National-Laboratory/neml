#include "interpolate.h"

#include "nemlmath.h"

#include <math.h>
#include <algorithm>
#include <limits>

namespace neml {

Interpolate::Interpolate() :
    valid_(true)
{

}

double Interpolate::operator()(double x) const
{
  return value(x);
}

bool Interpolate::valid() const
{
  return valid_;
}

InvalidInterpolate::InvalidInterpolate() :
    Interpolate()
{
  valid_ = false;
}

double InvalidInterpolate::value(double x) const
{
  return std::numeric_limits<double>::quiet_NaN();
}

PolynomialInterpolate::PolynomialInterpolate(const std::vector<double> coefs) :
    coefs_(coefs), Interpolate()
{

}

double PolynomialInterpolate::value(double x) const
{
  return polyval(&coefs_[0], coefs_.size(), x);
}


PiecewiseLinearInterpolate::PiecewiseLinearInterpolate(
    const std::vector<double> points,
    const std::vector<double> values) :
      points_(points), values_(values), Interpolate()
{
  // Check if sorted
  if (not std::is_sorted(points.begin(), points.end())) {
    valid_ = false; 
  }

  if (points.size() != values.size()) {
    valid_ = false;
  }
}

double PiecewiseLinearInterpolate::value(double x) const
{
  if (x < points_.front()) {
    return values_.front();
  }
  else if (x > points_.back()) {
    return values_.back();
  }
  else {
    auto it = points_.begin();
    for (it; it != points_.end(); ++it) {
      if (x <= *it) break;
    }
    size_t ind = std::distance(points_.begin(), it);
    double x1 = points_[ind-1];
    double x2 = points_[ind];
    double y1 = values_[ind-1];
    double y2 = values_[ind];

    return (y2-y1)/(x2-x1) * (x - x1) + y1;
  }
}

ConstantInterpolate::ConstantInterpolate(double v) :
    v_(v), Interpolate()
{

}

double ConstantInterpolate::value(double x) const
{
  return v_;
}

MTSShearInterpolate::MTSShearInterpolate(double V0, double D, double T0) :
    V0_(V0), D_(D), T0_(T0), Interpolate()
{
  
}

double MTSShearInterpolate::value(double x) const
{
  return V0_ - D_ / (exp(T0_ / x) - 1.0);
}

std::vector<std::shared_ptr<Interpolate>> 
    make_vector(const std::vector<double> & iv)
{
  std::vector<std::shared_ptr<Interpolate>> vt;
  for (auto it = iv.begin(); it != iv.end(); ++it) {
    vt.emplace_back(std::make_shared<ConstantInterpolate>(*it));
  }
  return vt;
}

std::vector<double> eval_vector(
    const std::vector<std::shared_ptr<Interpolate>> & iv, double x)
{
  std::vector<double> vt;
  for (auto it = iv.begin(); it != iv.end(); ++it) {
    vt.push_back((*it)->value(x));
  }
  return vt;
}

} // namespace neml
