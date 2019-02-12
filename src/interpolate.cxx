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

PolynomialInterpolate::PolynomialInterpolate(const std::vector<double> coefs) :
    Interpolate(), coefs_(coefs)
{
  int n = coefs_.size();
  deriv_.resize(n - 1);
  for (int i = 0; i < n - 1; i++) {
    deriv_[i] = coefs_[i] * ((double) (n - 1 - i));
  }
}

std::string PolynomialInterpolate::type()
{
  return "PolynomialInterpolate";
}

ParameterSet PolynomialInterpolate::parameters()
{
  ParameterSet pset(PolynomialInterpolate::type());

  pset.add_parameter<std::vector<double>>("coefs");

  return pset;
}

std::unique_ptr<NEMLObject> PolynomialInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<PolynomialInterpolate>(
      params.get_parameter<std::vector<double>>("coefs")
      ); 
}

double PolynomialInterpolate::value(double x) const
{
  return polyval(&coefs_[0], coefs_.size(), x);
}

double PolynomialInterpolate::derivative(double x) const
{
  return polyval(&deriv_[0], deriv_.size(), x);
}


PiecewiseLinearInterpolate::PiecewiseLinearInterpolate(
    const std::vector<double> points,
    const std::vector<double> values) :
      Interpolate(), points_(points), values_(values)
{
  // Check if sorted
  if (not std::is_sorted(points.begin(), points.end())) {
    valid_ = false; 
  }

  if (points.size() != values.size()) {
    valid_ = false;
  }
}

std::string PiecewiseLinearInterpolate::type()
{
  return "PiecewiseLinearInterpolate";
}

ParameterSet PiecewiseLinearInterpolate::parameters()
{
  ParameterSet pset(PiecewiseLinearInterpolate::type());

  pset.add_parameter<std::vector<double>>("points");
  pset.add_parameter<std::vector<double>>("values");

  return pset;
}

std::unique_ptr<NEMLObject> PiecewiseLinearInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<PiecewiseLinearInterpolate>(
      params.get_parameter<std::vector<double>>("points"),
      params.get_parameter<std::vector<double>>("values")
      ); 
}

double PiecewiseLinearInterpolate::value(double x) const
{
  if (x <= points_.front()) {
    return values_.front();
  }
  else if (x >= points_.back()) {
    return values_.back();
  }
  else {
    auto it = points_.begin();
    for (; it != points_.end(); ++it) {
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

double PiecewiseLinearInterpolate::derivative(double x) const
{
  if (x <= points_.front()) {
    return 0.0;
  }
  else if (x >= points_.back()) {
    return 0.0;
  }
  else {
    auto it = points_.begin();
    for (; it != points_.end(); ++it) {
      if (x <= *it) break;
    }
    size_t ind = std::distance(points_.begin(), it);
    double x1 = points_[ind-1];
    double x2 = points_[ind];
    double y1 = values_[ind-1];
    double y2 = values_[ind];

    return (y2-y1)/(x2-x1);
  }
}

PiecewiseLogLinearInterpolate::PiecewiseLogLinearInterpolate(
    const std::vector<double> points,
    const std::vector<double> values) :
      Interpolate(), points_(points), values_(values)
{
  // Check if sorted
  if (not std::is_sorted(points.begin(), points.end())) {
    valid_ = false; 
  }

  if (points.size() != values.size()) {
    valid_ = false;
  }

  for (auto it = values_.begin(); it != values_.end(); ++it) {
    if (*it < 0.0) valid_ = false;
    *it = log(*it);
  }
}

std::string PiecewiseLogLinearInterpolate::type()
{
  return "PiecewiseLogLinearInterpolate";
}

ParameterSet PiecewiseLogLinearInterpolate::parameters()
{
  ParameterSet pset(PiecewiseLogLinearInterpolate::type());

  pset.add_parameter<std::vector<double>>("points");
  pset.add_parameter<std::vector<double>>("values");

  return pset;
}

std::unique_ptr<NEMLObject> PiecewiseLogLinearInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<PiecewiseLogLinearInterpolate>(
      params.get_parameter<std::vector<double>>("points"),
      params.get_parameter<std::vector<double>>("values")
      ); 
}

double PiecewiseLogLinearInterpolate::value(double x) const
{
  if (x <= points_.front()) {
    return exp(values_.front());
  }
  else if (x >= points_.back()) {
    return exp(values_.back());
  }
  else {
    auto it = points_.begin();
    for (; it != points_.end(); ++it) {
      if (x <= *it) break;
    }
    size_t ind = std::distance(points_.begin(), it);
    double x1 = points_[ind-1];
    double x2 = points_[ind];
    double y1 = values_[ind-1];
    double y2 = values_[ind];

    return exp((y2-y1)/(x2-x1) * (x - x1) + y1);
  }
}

double PiecewiseLogLinearInterpolate::derivative(double x) const
{
  if (x <= points_.front()) {
    return 0.0;
  }
  else if (x >= points_.back()) {
    return 0.0;
  }
  else {
    auto it = points_.begin();
    for (; it != points_.end(); ++it) {
      if (x <= *it) break;
    }
    size_t ind = std::distance(points_.begin(), it);
    double x1 = points_[ind-1];
    double x2 = points_[ind];
    double y1 = values_[ind-1];
    double y2 = values_[ind];

    return exp((y2-y1)/(x2-x1) * (x-x1) + y1) * (y2-y1)/(x2-x1);
  }
}

ConstantInterpolate::ConstantInterpolate(double v) :
    Interpolate(), v_(v)
{

}

std::string ConstantInterpolate::type()
{
  return "ConstantInterpolate";
}

ParameterSet ConstantInterpolate::parameters()
{
  ParameterSet pset(ConstantInterpolate::type());

  pset.add_parameter<double>("v");

  return pset;
}

std::unique_ptr<NEMLObject> ConstantInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<ConstantInterpolate>(
      params.get_parameter<double>("v")
      ); 
}

double ConstantInterpolate::value(double x) const
{
  return v_;
}

double ConstantInterpolate::derivative(double x) const
{
  return 0.0;
}

MTSShearInterpolate::MTSShearInterpolate(double V0, double D, double T0) :
    Interpolate(), V0_(V0), D_(D), T0_(T0)
{
  
}

std::string MTSShearInterpolate::type()
{
  return "MTSShearInterpolate";
}

ParameterSet MTSShearInterpolate::parameters()
{
  ParameterSet pset(MTSShearInterpolate::type());

  pset.add_parameter<double>("V0");
  pset.add_parameter<double>("D");
  pset.add_parameter<double>("T0");

  return pset;
}

std::unique_ptr<NEMLObject> MTSShearInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<MTSShearInterpolate>(
      params.get_parameter<double>("V0"),
      params.get_parameter<double>("D"),
      params.get_parameter<double>("T0")
      ); 
}

double MTSShearInterpolate::value(double x) const
{
  return V0_ - D_ / (exp(T0_ / x) - 1.0);
}

double MTSShearInterpolate::derivative(double x) const
{
  return -D_ * T0_ / (4.0 * pow(x * sinh(T0_ / (2 * x)),2));
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

std::vector<double> eval_deriv_vector(
    const std::vector<std::shared_ptr<Interpolate>> & iv, double x)
{
  std::vector<double> vt;
  for (auto it = iv.begin(); it != iv.end(); ++it) {
    vt.push_back((*it)->derivative(x));
  }
  return vt;
}

} // namespace neml
