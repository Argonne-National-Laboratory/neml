#include "creep.h"

#include "math/nemlmath.h"
#include "nemlerror.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace neml {

ScalarCreepRule::ScalarCreepRule(ParameterSet & params) :
    NEMLObject(params)
{

}

// Scalar creep default derivatives for time and temperature
double ScalarCreepRule::dg_dt(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

double ScalarCreepRule::dg_dT(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

// Implementation of power law creep
PowerLawCreep::PowerLawCreep(ParameterSet & params) :
    ScalarCreepRule(params),
    A_(params.get_object_parameter<Interpolate>("A")), 
    n_(params.get_object_parameter<Interpolate>("n"))
{

}

std::string PowerLawCreep::type()
{
  return "PowerLawCreep";
}

ParameterSet PowerLawCreep::parameters()
{
  ParameterSet pset(PowerLawCreep::type());

  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

std::unique_ptr<NEMLObject> PowerLawCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<PowerLawCreep>(params);
}


double PowerLawCreep::g(double seq, double eeq, double t, double T) const
{
  return A_->value(T) * pow(seq, n_->value(T));
}

double PowerLawCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  double nv = n_->value(T);

  return A_->value(T) * nv * pow(seq, nv - 1.0);
}

double PowerLawCreep::dg_de(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

double PowerLawCreep::A(double T) const
{
  return A_->value(T);
}

double PowerLawCreep::n(double T) const
{
  return n_->value(T);
}

// Implementation of power law creep
NormalizedPowerLawCreep::NormalizedPowerLawCreep(ParameterSet & params) :
    ScalarCreepRule(params),
    s0_(params.get_object_parameter<Interpolate>("s0")), 
    n_(params.get_object_parameter<Interpolate>("n"))
{

}

std::string NormalizedPowerLawCreep::type()
{
  return "NormalizedPowerLawCreep";
}

ParameterSet NormalizedPowerLawCreep::parameters()
{
  ParameterSet pset(NormalizedPowerLawCreep::type());

  pset.add_parameter<NEMLObject>("s0");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

std::unique_ptr<NEMLObject> NormalizedPowerLawCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<NormalizedPowerLawCreep>(params); 
}


double NormalizedPowerLawCreep::g(double seq, double eeq, double t, double T) const
{
  return pow(seq / s0_->value(T), n_->value(T));
}

double NormalizedPowerLawCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  double nv = n_->value(T);
  double s0v = s0_->value(T);

  return nv / s0v * pow(seq / s0v, nv - 1.0);
}

double NormalizedPowerLawCreep::dg_de(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

// Implementation of the Blackburn minimum creep rate equation
BlackburnMinimumCreep::BlackburnMinimumCreep(ParameterSet & params) :
    ScalarCreepRule(params),
    A_(params.get_object_parameter<Interpolate>("A")),
    n_(params.get_object_parameter<Interpolate>("n")),
    beta_(params.get_object_parameter<Interpolate>("beta")),
    R_(params.get_parameter<double>("R")),
    Q_(params.get_parameter<double>("Q"))
{

}

std::string BlackburnMinimumCreep::type()
{
  return "BlackburnMinimumCreep";
}

ParameterSet BlackburnMinimumCreep::parameters()
{
  ParameterSet pset(BlackburnMinimumCreep::type());

  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("n");
  pset.add_parameter<NEMLObject>("beta");
  pset.add_parameter<double>("R");
  pset.add_parameter<double>("Q");

  return pset;
}

std::unique_ptr<NEMLObject> BlackburnMinimumCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<BlackburnMinimumCreep>(params); 
}


double BlackburnMinimumCreep::g(double seq, double eeq, double t, double T) const
{
  double A = A_->value(T);
  double n = n_->value(T);
  double beta = beta_->value(T);

  return A * pow(sinh(beta * seq / n), n) * exp(-Q_ / (R_ * T));
}

double BlackburnMinimumCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  double A = A_->value(T);
  double n = n_->value(T);
  double beta = beta_->value(T);

  return A * beta * exp(-Q_ / (R_ * T)) * cosh(beta * seq / n) * 
      pow(sinh(beta * seq / n), n - 1.0);

}

double BlackburnMinimumCreep::dg_de(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

double BlackburnMinimumCreep::dg_dT(double seq, double eeq, double t, double T) const
{
  double A = A_->value(T);
  double n = n_->value(T);
  double beta = beta_->value(T);

  return A * pow(sinh(beta * seq / n), n) * exp(-Q_ / (R_ * T)) * Q_ / (R_ * T * T);
}

// Implementation of the Swindeman minimum creep rate equation
SwindemanMinimumCreep::SwindemanMinimumCreep(ParameterSet & params) :
    ScalarCreepRule(params),
    C_(params.get_parameter<double>("C")), 
    n_(params.get_parameter<double>("n")),
    V_(params.get_parameter<double>("V")),
    Q_(params.get_parameter<double>("Q")),
    shift_(params.get_parameter<bool>("celsius") ? 273.15 : 0.0)
{

}

std::string SwindemanMinimumCreep::type()
{
  return "SwindemanMinimumCreep";
}

ParameterSet SwindemanMinimumCreep::parameters()
{
  ParameterSet pset(SwindemanMinimumCreep::type());

  pset.add_parameter<double>("C");
  pset.add_parameter<double>("n");
  pset.add_parameter<double>("V");
  pset.add_parameter<double>("Q");

  pset.add_optional_parameter<bool>("celsius", false);

  return pset;
}

std::unique_ptr<NEMLObject> SwindemanMinimumCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<SwindemanMinimumCreep>(params); 
}


double SwindemanMinimumCreep::g(double seq, double eeq, double t, double T) const
{
  return C_ * pow(seq, n_) * exp(V_ * seq) * exp(-Q_/(T+shift_));
}

double SwindemanMinimumCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  return C_ * exp(-Q_/(T+shift_)) * (n_ + seq * V_) * exp(seq * V_) * pow(seq, n_ - 1.0);
}

double SwindemanMinimumCreep::dg_de(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

double SwindemanMinimumCreep::dg_dT(double seq, double eeq, double t, double T) const
{
  return C_ * pow(seq, n_) * exp(V_ * seq) * exp(-Q_/(T+shift_)) * Q_ /
      ((T+shift_) * (T+shift_)); 
}

// Implementation of the mechanism switching model
RegionKMCreep::RegionKMCreep(ParameterSet & params) :
    ScalarCreepRule(params),
    cuts_(params.get_parameter<std::vector<double>>("cuts")),
    A_(params.get_object_parameter_vector<Interpolate>("A")), 
    B_(params.get_object_parameter_vector<Interpolate>("B")),
    kboltz_(params.get_parameter<double>("kboltz")),
    b_(params.get_parameter<double>("b")),
    eps0_( params.get_parameter<double>("eps0")), 
    b3_(pow(b_,3)),
    emodel_(params.get_object_parameter<LinearElasticModel>("emodel")),
    shift_(params.get_parameter<bool>("celsius") ? 273.15 : 0.0)
{

}

std::string RegionKMCreep::type()
{
  return "RegionKMCreep";
}

ParameterSet RegionKMCreep::parameters()
{
  ParameterSet pset(RegionKMCreep::type());

  pset.add_parameter<std::vector<double>>("cuts");
  pset.add_parameter<std::vector<NEMLObject>>("A");
  pset.add_parameter<std::vector<NEMLObject>>("B");
  pset.add_parameter<double>("kboltz");
  pset.add_parameter<double>("b");
  pset.add_parameter<double>("eps0");
  pset.add_parameter<NEMLObject>("emodel");

  pset.add_optional_parameter<bool>("celsius", false);

  return pset;
}

std::unique_ptr<NEMLObject> RegionKMCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<RegionKMCreep>(params); 
}

double RegionKMCreep::g(double seq, double eeq, double t, double T) const
{
  double A, B;
  select_region_(seq, T, A, B);
  double G = emodel_->G(T);
  double C1 = -G * b3_ / (kboltz_*(T+shift_));

  return eps0_ * exp(C1*B) * pow(seq / G, C1 * A);
}

double RegionKMCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  double A, B;
  select_region_(seq, T, A, B);
  double G = emodel_->G(T);
  double C1 = -G * b3_ / (kboltz_*(T+shift_));
  
  return eps0_ * exp(C1*B) * C1 * A / G * pow(seq / G, C1 * A - 1.0);
}

double RegionKMCreep::dg_de(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

void RegionKMCreep::select_region_(double seq, double T, double & Ai, double & Bi) const
{
  double mu = emodel_->G(T);
  double neq = seq / mu;
  int nregion = A_.size();
  if (nregion == 1) {
    Ai = A_[0]->value(T);
    Bi = B_[0]->value(T);
    return;
  }

  size_t i;
  if (neq < cuts_[0]) {
    Ai = A_[0]->value(T);
    Bi = B_[0]->value(T);
    return;
  }
  for (i=0; i<cuts_.size(); i++) {
    if (neq > cuts_[i]) {
      Ai = A_[i+1]->value(T);
      Bi = B_[i+1]->value(T);
      return;
    }
  }
  if (i == cuts_.size()) {
    Ai = A_[i]->value(T);
    Bi = B_[i]->value(T);
  }
}

// Implementation of Norton-Bailey creep
NortonBaileyCreep::NortonBaileyCreep(ParameterSet & params) :
    ScalarCreepRule(params),
    A_(params.get_object_parameter<Interpolate>("A")),
    m_(params.get_object_parameter<Interpolate>("m")),
    n_(params.get_object_parameter<Interpolate>("n"))
{

}

std::string NortonBaileyCreep::type()
{
  return "NortonBaileyCreep";
}

ParameterSet NortonBaileyCreep::parameters()
{
  ParameterSet pset(NortonBaileyCreep::type());

  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("m");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

std::unique_ptr<NEMLObject> NortonBaileyCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<NortonBaileyCreep>(params); 
}

double NortonBaileyCreep::g(double seq, double eeq, double t, double T) const
{
  double A = A_->value(T);
  double m = m_->value(T);
  double n = n_->value(T);

  // Hack, really should figure out limits
  if (seq < std::numeric_limits<double>::epsilon()) {
    seq = std::numeric_limits<double>::epsilon(); 
  }
  if (eeq < std::numeric_limits<double>::epsilon()) {
    eeq = std::numeric_limits<double>::epsilon(); 
  }

  return m * pow(A, 1.0 / m) * pow(seq, n / m) * pow(eeq, (m - 1.0) / m); 
}

double NortonBaileyCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  double A = A_->value(T);
  double m = m_->value(T);
  double n = n_->value(T);

  // Hack, really should figure out limits
  if (seq < std::numeric_limits<double>::epsilon()) {
    seq = std::numeric_limits<double>::epsilon(); 
  }
  if (eeq < std::numeric_limits<double>::epsilon()) {
    eeq = std::numeric_limits<double>::epsilon(); 
  }

  return n * pow(A, 1.0 / m) * pow(seq, n / m - 1.0) * pow(eeq, (m - 1.0) / m);
}

double NortonBaileyCreep::dg_de(double seq, double eeq, double t, double T) const
{
  double A = A_->value(T);
  double m = m_->value(T);
  double n = n_->value(T);

  // Hack, really should figure out limits
  if (seq < std::numeric_limits<double>::epsilon()) {
    seq = std::numeric_limits<double>::epsilon(); 
  }
  if (eeq < std::numeric_limits<double>::epsilon()) {
    eeq = std::numeric_limits<double>::epsilon(); 
  }

  return (m - 1) * pow(A, 1.0 / m) * pow(seq, n / m) * pow(eeq, -1.0 / m);
}

double NortonBaileyCreep::A(double T) const
{
  return A_->value(T);
}

double NortonBaileyCreep::m(double T) const
{
  return m_->value(T);
}

double NortonBaileyCreep::n(double T) const
{
  return n_->value(T);
}


MukherjeeCreep::MukherjeeCreep(ParameterSet & params) :
    ScalarCreepRule(params),
    emodel_(params.get_object_parameter<LinearElasticModel>("emodel")), 
    A_(params.get_parameter<double>("A")),
    n_(params.get_parameter<double>("n")),
    D0_(params.get_parameter<double>("D0")), 
    Q_(params.get_parameter<double>("Q")),
    b_(params.get_parameter<double>("b")),
    k_(params.get_parameter<double>("k")), 
    R_(params.get_parameter<double>("R"))
{
  
}

std::string MukherjeeCreep::type()
{
  return "MukherjeeCreep";
}

ParameterSet MukherjeeCreep::parameters()
{
  ParameterSet pset(MukherjeeCreep::type());

  pset.add_parameter<NEMLObject>("emodel");
  pset.add_parameter<double>("A");
  pset.add_parameter<double>("n");
  pset.add_parameter<double>("D0");
  pset.add_parameter<double>("Q");
  pset.add_parameter<double>("b");
  pset.add_parameter<double>("k");
  pset.add_parameter<double>("R");

  return pset;
}

std::unique_ptr<NEMLObject> MukherjeeCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<MukherjeeCreep>(params); 
}

double MukherjeeCreep::g(double seq, double eeq, double t, double T) const
{
  double mu = emodel_->G(T);
  double Dv = D0_ * exp(-Q_ / (R_ * T));
  return A_ * Dv * mu * b_ / (k_ * T) * pow(seq / mu, n_);
}

double MukherjeeCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  double mu = emodel_->G(T);
  double Dv = D0_ * exp(-Q_ / (R_ * T));
  return n_ * A_ * Dv * mu * b_ / (k_ * T) * pow(seq / mu, n_ - 1.0) / mu;
}

double MukherjeeCreep::dg_de(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

double MukherjeeCreep::A() const
{
  return A_;
}

double MukherjeeCreep::n() const
{
  return n_;
}

double MukherjeeCreep::D0() const
{
  return D0_;
}

double MukherjeeCreep::Q() const
{
  return Q_;
}

double MukherjeeCreep::b() const
{
  return b_;
}

double MukherjeeCreep::k() const
{
  return k_;
}

double MukherjeeCreep::R() const
{
  return R_;
}

GenericCreep::GenericCreep(ParameterSet & params) :
    ScalarCreepRule(params),
    cfn_(params.get_object_parameter<Interpolate>("cfn"))
{

}

std::string GenericCreep::type()
{
  return "GenericCreep";
}

std::unique_ptr<NEMLObject> GenericCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<GenericCreep>(params);
}

ParameterSet GenericCreep::parameters()
{
  ParameterSet pset(GenericCreep::type());

  pset.add_parameter<NEMLObject>("cfn");

  return pset;
}

double GenericCreep::g(double seq, double eeq, double t, double T) const
{
  return exp(cfn_->value(log(seq)));
}

double GenericCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  double f = cfn_->value(log(seq));
  double df = cfn_->derivative(log(seq));
  
  if (seq > 0.0) {
    return f * exp(f) * df / seq;
  }
  else {
    return 0.0;
  }
}

double GenericCreep::dg_de(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

// Setup for solve
CreepModel::CreepModel(ParameterSet & params) :
    NEMLObject(params),
    rtol_(params.get_parameter<double>("rtol")),
    atol_(params.get_parameter<double>("atol")),
    miter_(params.get_parameter<int>("miter")),
    verbose_(params.get_parameter<bool>("verbose")), 
    linesearch_(params.get_parameter<bool>("linesearch"))
{

}

// Creep model default derivatives for time and temperature
Symmetric CreepModel::df_dt(const Symmetric & s, const Symmetric & e, double t, double T) const
{
  return Symmetric();
}

Symmetric CreepModel::df_dT(const Symmetric & s, const Symmetric & e, double t, double T) const
{
  return Symmetric();
}

// Implementation of creep model update
void CreepModel::update(const Symmetric & s_np1, 
                        Symmetric & e_np1, const Symmetric & e_n,
                        double T_np1, double T_n,
                        double t_np1, double t_n,
                        SymSymR4 & A_np1)
{
  // Setup the trial state
  std::unique_ptr<CreepModelTrialState> ts = make_trial_state(s_np1, e_n, T_np1, T_n, t_np1, t_n);

  // Solve for the new creep strain
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  solve(this, x, ts.get(), {rtol_, atol_, miter_, verbose_, linesearch_});
  
  // Extract
  std::copy(x, x+6, e_np1.s());

  // Get the tangent
  A_np1 = calc_tangent_(e_np1, *ts);
}

std::unique_ptr<CreepModelTrialState> CreepModel::make_trial_state(
    const Symmetric & s_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n) const
{
  return std::make_unique<CreepModelTrialState>(s_np1, e_n, T_np1, t_np1 - t_n, t_np1);
}

// Implement the solve
size_t CreepModel::nparams() const
{
  return 6; // the creep strain
}

void CreepModel::init_x(double * const x, TrialState * ts)
{
  CreepModelTrialState * tss = static_cast<CreepModelTrialState *>(ts);

  // Just make it the previous value
  std::copy(tss->e_n.data(), tss->e_n.data()+6, x);
}

void CreepModel::RJ(const double * const x, TrialState * ts, 
                     double * const R, double * const J)
{
  CreepModelTrialState * tss = static_cast<CreepModelTrialState *>(ts);
  
  Symmetric e_np1(x);
  Symmetric R1(R);
  SymSymR4 J1(J);

  // Residual
  R1 = e_np1 - tss->e_n - f(tss->s_np1, e_np1, tss->t, tss->T) * tss->dt;

  // Jacobian
  J1 = SymSymR4::id() - df_de(tss->s_np1, e_np1, tss->t, tss->T) *tss->dt;
}

// Helper for tangent
SymSymR4 CreepModel::calc_tangent_(const Symmetric & e_np1, CreepModelTrialState & ts)
{
  Symmetric R;
  SymSymR4 J;

  RJ(e_np1.data(), &ts, R.s(), J.s());

  SymSymR4 B = df_ds(ts.s_np1, e_np1, ts.t, ts.T);

  return (J.inverse() * ts.dt).dot(B);
}

// Implementation of 2.25Cr-1Mo rule
MinCreep225Cr1MoCreep::MinCreep225Cr1MoCreep(ParameterSet & params) :
    ScalarCreepRule(params)
{

}

std::string MinCreep225Cr1MoCreep::type()
{
  return "MinCreep225Cr1MoCreep";
}

ParameterSet MinCreep225Cr1MoCreep::parameters()
{
  ParameterSet pset(MinCreep225Cr1MoCreep::type());

  return pset;
}

std::unique_ptr<NEMLObject> MinCreep225Cr1MoCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<MinCreep225Cr1MoCreep>(params); 
}


double MinCreep225Cr1MoCreep::g(double seq, double eeq, double t, double T) const
{
  if (seq < 60.0) {
    return e1_(seq, T);
  }
  else {
    if (T <= (13.571 * pow(seq, 0.68127) - 1.8 * seq + 710.78)) {
      return e1_(seq, T);
    }
    else {
      return e2_(seq, T);
    }
  }

}

double MinCreep225Cr1MoCreep::dg_ds(double seq, double eeq, double t, double T) const
{
  if (seq < 60.0) {
    return de1_(seq, T);
  }
  else {
    if (T <= (13.571 * pow(seq, 0.68127) - 1.8 * seq + 710.78)) {
      return de1_(seq, T);
    }
    else {
      return de2_(seq, T);
    }
  }

}

double MinCreep225Cr1MoCreep::dg_de(double seq, double eeq, double t, double T) const
{
  return 0.0;
}

double MinCreep225Cr1MoCreep::e1_(double seq, double T) const
{
  double U = MinCreep225Cr1MoCreep::U->value(T);
  double exp = 6.7475 + 0.011426 * seq + 987.72 / U * log10(seq) - 13494.0/T;
  return pow(10.0, exp) / 100.0;
}

double MinCreep225Cr1MoCreep::e2_(double seq, double T) const
{
  double U = MinCreep225Cr1MoCreep::U->value(T);
  double exp = 11.498 - 8.2226*U / T - 20448 / T + 5862.4 / T * log10(seq);
  return pow(10.0, exp) / 100.0;
}

double MinCreep225Cr1MoCreep::de1_(double seq, double T) const
{
  double U = MinCreep225Cr1MoCreep::U->value(T);
  double exp = 4.7475 + 0.011426 * seq + 428.961 / U * log(seq) - 13494.0/T;
  return pow(10.0, exp) * (0.011426 + 428.961 / (seq * U)) * log(10.0);
}

double MinCreep225Cr1MoCreep::de2_(double seq, double T) const
{
  double U = MinCreep225Cr1MoCreep::U->value(T);
  double exp = 9.498 - 8.2226*U / T - 20448 / T + 2546.01 / T * log(seq);
  return 2546.01 * pow(10.0, exp) * log(10.0) / (seq * T);
}

const std::shared_ptr<PiecewiseLinearInterpolate> MinCreep225Cr1MoCreep::U = make_piecewise({644.15,673.15,723.15,773.15,823.15,873.15,894.15,922.15},{471,468,452,418,634,284,300,270});

// Implementation of J2 creep
J2CreepModel::J2CreepModel(ParameterSet & params) :
    CreepModel(params), rule_(params.get_object_parameter<ScalarCreepRule>("rule"))
{

}

std::string J2CreepModel::type()
{
  return "J2CreepModel";
}

ParameterSet J2CreepModel::parameters()
{
  ParameterSet pset(J2CreepModel::type());

  pset.add_parameter<NEMLObject>("rule");
  
  pset.add_optional_parameter<double>("rtol", 1.0e-8);
  pset.add_optional_parameter<double>("atol", 1.0e-10);
  pset.add_optional_parameter<int>("miter", 25);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<bool>("linesearch", false);

  return pset;
}

std::unique_ptr<NEMLObject> J2CreepModel::initialize(ParameterSet & params)
{
  return neml::make_unique<J2CreepModel>(params); 
}

Symmetric J2CreepModel::f(const Symmetric & s, const Symmetric & e, double t, double T) const
{
  // Rate * direction
  return 3.0/2.0 * sdir(s) * rule_->g(seq(s), eeq(e), t, T);
}

SymSymR4 J2CreepModel::df_ds(const Symmetric & s, const Symmetric & e, double t, double T) const
{
  // Gather the effective stresses
  double se = seq(s);
  double ee = eeq(e);
  
  // Get the direction
  Symmetric dir = sdir(s);

  // Get the rate
  double rate = rule_->g(se, ee, t, T);

  // Get the rate derivative
  double drate = rule_->dg_ds(se, ee, t, T);

  SymSymR4 ID = SymSymR4::id_dev() * 3.0 / 2.0;

  // Begin forming outer products
  // Hack, really should figure out limit se -> 0.0
  if (se < std::numeric_limits<double>::epsilon()) {
    se = std::numeric_limits<double>::epsilon(); 
  }

  SymSymR4 A = 3.0/2.0 * douter(dir, dir) * (drate - rate / se) + SymSymR4::id() * rate / se;
  
  return A.dot(ID);
}

SymSymR4 J2CreepModel::df_de(const Symmetric & s, const Symmetric & e, double t, double T) const
{
  // Get the derivative
  double drate = rule_->dg_de(seq(s), eeq(e), t, T);

  // Get the stress and strain direction
  Symmetric s_dir = sdir(s);
  Symmetric e_dir = edir(e) * drate;

  return douter(s_dir, e_dir);
}

Symmetric J2CreepModel::df_dt(const Symmetric & s, const Symmetric & e, double t, double T) const 
{
  return 3.0/2.0 * sdir(s) * rule_->dg_dt(seq(s), eeq(e), t, T);
}

Symmetric J2CreepModel::df_dT(const Symmetric & s, const Symmetric & e, double t, double T) const
{
  return 3.0/2.0 * sdir(s) * rule_->dg_dT(seq(s), eeq(e), t, T); 
}

// Helpers for J2 plasticity
double J2CreepModel::seq(const Symmetric & s) const
{
  return sqrt(3.0/2.0)*s.dev().norm();
}

double J2CreepModel::eeq(const Symmetric & e) const
{
  return sqrt(2.0 / 3.0) * e.norm();
}

Symmetric J2CreepModel::sdir(const Symmetric & s) const
{
  double se = seq(s);
  if (se < std::numeric_limits<double>::epsilon()) {
    return Symmetric();
  }
  else { 
    return s.dev() / se;
  }
}

Symmetric J2CreepModel::edir(const Symmetric & e) const
{
  double ee = eeq(e);
  if (ee < std::numeric_limits<double>::epsilon()) {
    return Symmetric();
  }
  else {
    return e / ee;
  }
}

} // namespace neml
