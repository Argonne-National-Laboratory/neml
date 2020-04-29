#include "interpolate.h"

#include "math/nemlmath.h"

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
  deriv_ = differentiate_poly(coefs_);
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
  return polyval(coefs_, x);
}

double PolynomialInterpolate::derivative(double x) const
{
  return polyval(deriv_, x);
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

GenericPiecewiseInterpolate::GenericPiecewiseInterpolate(
    std::vector<double> points,
    std::vector<std::shared_ptr<Interpolate>> functions) :
      Interpolate(), points_(points), functions_(functions)
{
  // Check if sorted
  if (not std::is_sorted(points.begin(), points.end())) {
    valid_ = false; 
  }

  if (points.size() != (functions.size()+1)) {
    valid_ = false;
  }
}

std::string GenericPiecewiseInterpolate::type()
{
  return "GenericPiecewiseInterpolate";
}

ParameterSet GenericPiecewiseInterpolate::parameters()
{
  ParameterSet pset(GenericPiecewiseInterpolate::type());

  pset.add_parameter<std::vector<double>>("points");
  pset.add_parameter<std::vector<NEMLObject>>("functions");

  return pset;
}

std::unique_ptr<NEMLObject> GenericPiecewiseInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<GenericPiecewiseInterpolate>(
      params.get_parameter<std::vector<double>>("points"),
      params.get_object_parameter_vector<Interpolate>("functions")
      ); 
}

double GenericPiecewiseInterpolate::value(double x) const
{
  if (x <= points_.front()) {
    return functions_[0]->value(x);
  }
  else if (x >= points_.back()) {
    return functions_.back()->value(x);
  }
  else {
    auto it = points_.begin();
    for (; it != points_.end(); ++it) {
      if (x <= *it) break;
    }
    size_t ind = std::distance(points_.begin(), it);

    return functions_[ind]->value(x);
  }
}

double GenericPiecewiseInterpolate::derivative(double x) const
{
  if (x <= points_.front()) {
    return functions_[0]->derivative(x);
  }
  else if (x >= points_.back()) {
    return functions_.back()->derivative(x);
  }
  else {
    auto it = points_.begin();
    for (; it != points_.end(); ++it) {
      if (x <= *it) break;
    }
    size_t ind = std::distance(points_.begin(), it);

    return functions_[ind]->derivative(x);
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

PiecewiseSemiLogLinearInterpolate::PiecewiseSemiLogLinearInterpolate(
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

  for (auto pt : points) {
    if (pt < 0.0) {
      valid_ = false;
    }
  }
}

std::string PiecewiseSemiLogLinearInterpolate::type()
{
  return "PiecewiseSemiLogLinearInterpolate";
}

ParameterSet PiecewiseSemiLogLinearInterpolate::parameters()
{
  ParameterSet pset(PiecewiseSemiLogLinearInterpolate::type());

  pset.add_parameter<std::vector<double>>("points");
  pset.add_parameter<std::vector<double>>("values");

  return pset;
}

std::unique_ptr<NEMLObject> PiecewiseSemiLogLinearInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<PiecewiseSemiLogLinearInterpolate>(
      params.get_parameter<std::vector<double>>("points"),
      params.get_parameter<std::vector<double>>("values")
      ); 
}

double PiecewiseSemiLogLinearInterpolate::value(double x) const
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

    return (y2-y1)/(log10(x2)-log10(x1)) * (log10(x) - log10(x1)) + y1;
  }
}

double PiecewiseSemiLogLinearInterpolate::derivative(double x) const
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

    return (y2-y1)/(log10(x2)-log10(x1)) / (x * log(10));
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

ExpInterpolate::ExpInterpolate(double A, double B) :
    Interpolate(), A_(A), B_(B)
{

}

std::string ExpInterpolate::type()
{
  return "ExpInterpolate";
}

ParameterSet ExpInterpolate::parameters()
{
  ParameterSet pset(ExpInterpolate::type());

  pset.add_parameter<double>("A");
  pset.add_parameter<double>("B");

  return pset;
}

std::unique_ptr<NEMLObject> ExpInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<ExpInterpolate>(
      params.get_parameter<double>("A"),
      params.get_parameter<double>("B")
      ); 
}

double ExpInterpolate::value(double x) const
{
  return A_*exp(B_/x);
}

double ExpInterpolate::derivative(double x) const
{
  return -A_ * B_ * exp(B_ / x) / (x*x);
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
