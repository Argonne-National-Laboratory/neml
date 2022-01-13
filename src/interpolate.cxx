#include "interpolate.h"

#include "math/nemlmath.h"

#include <math.h>
#include <algorithm>
#include <limits>

namespace neml {

Interpolate::Interpolate(ParameterSet & params) :
    NEMLObject(params),
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

PolynomialInterpolate::PolynomialInterpolate(ParameterSet & params) :
    Interpolate(params), 
    coefs_(params.get_parameter<std::vector<double>>("coefs"))
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
  return neml::make_unique<PolynomialInterpolate>(params); 
}

double PolynomialInterpolate::value(double x) const
{
  return polyval(coefs_, x);
}

double PolynomialInterpolate::derivative(double x) const
{
  return polyval(deriv_, x);
}


PiecewiseLinearInterpolate::PiecewiseLinearInterpolate(ParameterSet & params) :
      Interpolate(params),
      points_(params.get_parameter<std::vector<double>>("points")), 
      values_(params.get_parameter<std::vector<double>>("values"))
{
  // Check if sorted
  if (not std::is_sorted(points_.begin(), points_.end())) {
    valid_ = false; 
  }

  if (points_.size() != values_.size()) {
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
  return neml::make_unique<PiecewiseLinearInterpolate>(params); 
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

GenericPiecewiseInterpolate::GenericPiecewiseInterpolate(ParameterSet & params) :
      Interpolate(params),
      points_(params.get_parameter<std::vector<double>>("points")), 
      functions_(params.get_object_parameter_vector<Interpolate>("functions"))
{
  // Check if sorted
  if (not std::is_sorted(points_.begin(), points_.end())) {
    valid_ = false; 
  }

  if (points_.size() != (functions_.size()+1)) {
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
  return neml::make_unique<GenericPiecewiseInterpolate>(params); 
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
    ParameterSet & params) :
      Interpolate(params),
      points_(params.get_parameter<std::vector<double>>("points")),
      values_(params.get_parameter<std::vector<double>>("values"))
{
  // Check if sorted
  if (not std::is_sorted(points_.begin(), points_.end())) {
    valid_ = false; 
  }

  if (points_.size() != values_.size()) {
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
  return neml::make_unique<PiecewiseLogLinearInterpolate>(params); 
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

PiecewiseSemiLogXLinearInterpolate::PiecewiseSemiLogXLinearInterpolate(
    ParameterSet & params) :
      Interpolate(params), 
      points_(params.get_parameter<std::vector<double>>("points")), 
      values_(params.get_parameter<std::vector<double>>("values"))
{
  // Check if sorted
  if (not std::is_sorted(points_.begin(), points_.end())) {
    valid_ = false; 
  }

  if (points_.size() != values_.size()) {
    valid_ = false;
  }

  for (auto pt : points_) {
    if (pt < 0.0) {
      valid_ = false;
    }
  }
}

std::string PiecewiseSemiLogXLinearInterpolate::type()
{
  return "PiecewiseSemiLogXLinearInterpolate";
}

ParameterSet PiecewiseSemiLogXLinearInterpolate::parameters()
{
  ParameterSet pset(PiecewiseSemiLogXLinearInterpolate::type());

  pset.add_parameter<std::vector<double>>("points");
  pset.add_parameter<std::vector<double>>("values");

  return pset;
}

std::unique_ptr<NEMLObject> PiecewiseSemiLogXLinearInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<PiecewiseSemiLogXLinearInterpolate>(params); 
}

double PiecewiseSemiLogXLinearInterpolate::value(double x) const
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

double PiecewiseSemiLogXLinearInterpolate::derivative(double x) const
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

    return (y2-y1)/(std::log10(x2)-std::log10(x1)) / (x * std::log(10));
  }
}

ConstantInterpolate::ConstantInterpolate(ParameterSet & params) :
    Interpolate(params),
    v_(params.get_parameter<double>("v"))
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
  return neml::make_unique<ConstantInterpolate>(params); 
}

double ConstantInterpolate::value(double x) const
{
  return v_;
}

double ConstantInterpolate::derivative(double x) const
{
  return 0.0;
}

ExpInterpolate::ExpInterpolate(ParameterSet & params) :
    Interpolate(params),
    A_(params.get_parameter<double>("A")),
    B_(params.get_parameter<double>("B"))
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
  return neml::make_unique<ExpInterpolate>(params); 
}

double ExpInterpolate::value(double x) const
{
  return A_*exp(B_/x);
}

double ExpInterpolate::derivative(double x) const
{
  return -A_ * B_ * exp(B_ / x) / (x*x);
}

PowerLawInterpolate::PowerLawInterpolate(ParameterSet & params) :
    Interpolate(params),
    A_(params.get_parameter<double>("A")),
    B_(params.get_parameter<double>("B"))
{

}

std::string PowerLawInterpolate::type()
{
  return "PowerLawInterpolate";
}

ParameterSet PowerLawInterpolate::parameters()
{
  ParameterSet pset(PowerLawInterpolate::type());

  pset.add_parameter<double>("A");
  pset.add_parameter<double>("B");

  return pset;
}

std::unique_ptr<NEMLObject> PowerLawInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<PowerLawInterpolate>(params); 
}

double PowerLawInterpolate::value(double x) const
{
  return A_ * std::pow(x, B_);
}

double PowerLawInterpolate::derivative(double x) const
{
  return A_ * B_ * std::pow(x, B_-1.0);
}

MTSShearInterpolate::MTSShearInterpolate(ParameterSet & params) :
    Interpolate(params),
    V0_(params.get_parameter<double>("V0")), 
    D_(params.get_parameter<double>("D")), 
    T0_(params.get_parameter<double>("T0"))
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
  return neml::make_unique<MTSShearInterpolate>(params); 
}

double MTSShearInterpolate::value(double x) const
{
  return V0_ - D_ / (exp(T0_ / x) - 1.0);
}

double MTSShearInterpolate::derivative(double x) const
{
  return -D_ * T0_ / (4.0 * pow(x * sinh(T0_ / (2 * x)),2));
}


MTSInterpolate::MTSInterpolate(ParameterSet & params) :
    Interpolate(params),
    tau0_(params.get_parameter<double>("tau0")), 
    g0_(params.get_parameter<double>("g0")), 
    q_(params.get_parameter<double>("q")),
    p_(params.get_parameter<double>("p")),
    k_(params.get_parameter<double>("k")),
    b_(params.get_parameter<double>("b")),
    mu_(params.get_object_parameter<Interpolate>("mu"))
{
  
}

std::string MTSInterpolate::type()
{
  return "MTSInterpolate";
}

ParameterSet MTSInterpolate::parameters()
{
  ParameterSet pset(MTSInterpolate::type());

  pset.add_parameter<double>("tau0");
  pset.add_parameter<double>("g0");
  pset.add_parameter<double>("q");
  pset.add_parameter<double>("p");
  pset.add_parameter<double>("k");
  pset.add_parameter<double>("b");
  pset.add_parameter<NEMLObject>("mu");

  return pset;
}

std::unique_ptr<NEMLObject> MTSInterpolate::initialize(ParameterSet & params)
{
  return neml::make_unique<MTSInterpolate>(params); 
}

double MTSInterpolate::value(double x) const
{
  return tau0_ * std::pow(1.0 - 
                          std::pow(k_*x/(mu_->value(x) * std::pow(b_,3) * g0_),
                                   1.0/q_)
                          , 1.0/p_);
}

double MTSInterpolate::derivative(double x) const
{
  double b3 = std::pow(b_,3);
  double mu = mu_->value(x);
  double dmu = mu_->derivative(x);
  double inner = k_*x/(b3*g0_*mu);

  double A = -tau0_ * std::pow(inner, 1/q_) * std::pow(1-std::pow(inner, 1/q_),
                                                       1/p_ - 1) / (p_ * q_ * x);
  double B =  tau0_ * std::pow(inner, 1/q_) * std::pow(1-std::pow(inner, 1/q_), 
                                                       1/p_ - 1) / (p_ * q_ *
                                                                   mu);
  return A + B*dmu;
}

std::vector<std::shared_ptr<Interpolate>> 
    make_vector(const std::vector<double> & iv)
{
  std::vector<std::shared_ptr<Interpolate>> vt;
  for (auto it = iv.begin(); it != iv.end(); ++it) {
    vt.emplace_back(make_constant(*it));
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

std::shared_ptr<ConstantInterpolate> make_constant(double v)
{
  ParameterSet params = ConstantInterpolate::parameters();
  params.assign_parameter("v", v);
  
  return std::make_shared<ConstantInterpolate>(params);
}

std::unique_ptr<ConstantInterpolate> make_constant_unique(double v)
{
  ParameterSet params = ConstantInterpolate::parameters();
  params.assign_parameter("v", v);

  return neml::make_unique<ConstantInterpolate>(params);
}

std::shared_ptr<PiecewiseLinearInterpolate> make_piecewise(
    std::vector<double> points, std::vector<double> values)
{
  ParameterSet params = PiecewiseLinearInterpolate::parameters();
  params.assign_parameter("points", points);
  params.assign_parameter("values", values);
  
  return std::make_shared<PiecewiseLinearInterpolate>(params);
}



} // namespace neml
