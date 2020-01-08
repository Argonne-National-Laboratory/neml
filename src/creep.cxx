#include "creep.h"

#include "math/nemlmath.h"
#include "nemlerror.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace neml {

// Scalar creep default derivatives for time and temperature
int ScalarCreepRule::dg_dt(double seq, double eeq, double t, double T,
                           double & dg) const
{
  dg = 0.0;

  return 0;
}

int ScalarCreepRule::dg_dT(double seq, double eeq, double t, double T,
                           double & dg) const
{
  dg = 0.0;

  return 0;
}

// Implementation of power law creep
PowerLawCreep::PowerLawCreep(std::shared_ptr<Interpolate> A,
                             std::shared_ptr<Interpolate> n) :
    A_(A), n_(n)
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
  return neml::make_unique<PowerLawCreep>(
      params.get_object_parameter<Interpolate>("A"),
      params.get_object_parameter<Interpolate>("n")
      ); 
}


int PowerLawCreep::g(double seq, double eeq, double t, double T, double & g) const
{
  g = A_->value(T) * pow(seq, n_->value(T));
  return 0;
}

int PowerLawCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
{
  double nv = n_->value(T);

  dg = A_->value(T) * nv * pow(seq, nv - 1.0);
  return 0;
}

int PowerLawCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
{
  dg = 0.0;
  return 0;
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
NormalizedPowerLawCreep::NormalizedPowerLawCreep(std::shared_ptr<Interpolate> s0,
                             std::shared_ptr<Interpolate> n) :
    s0_(s0), n_(n)
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
  return neml::make_unique<NormalizedPowerLawCreep>(
      params.get_object_parameter<Interpolate>("s0"),
      params.get_object_parameter<Interpolate>("n")
      ); 
}


int NormalizedPowerLawCreep::g(double seq, double eeq, double t, double T, double & g) const
{
  g = pow(seq / s0_->value(T), n_->value(T));
  return 0;
}

int NormalizedPowerLawCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
{
  double nv = n_->value(T);
  double s0v = s0_->value(T);

  dg = nv / s0v * pow(seq / s0v, nv - 1.0);
  return 0;
}

int NormalizedPowerLawCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
{
  dg = 0.0;
  return 0;
}

// Implementation of the Blackburn minimum creep rate equation
BlackburnMinimumCreep::BlackburnMinimumCreep(
    std::shared_ptr<Interpolate> A,
    std::shared_ptr<Interpolate> n,
    std::shared_ptr<Interpolate> beta,
    double R, double Q) :
      A_(A), n_(n), beta_(beta), R_(R), Q_(Q)
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
  return neml::make_unique<BlackburnMinimumCreep>(
      params.get_object_parameter<Interpolate>("A"),
      params.get_object_parameter<Interpolate>("n"),
      params.get_object_parameter<Interpolate>("beta"),
      params.get_parameter<double>("R"),
      params.get_parameter<double>("Q")
      ); 
}


int BlackburnMinimumCreep::g(double seq, double eeq, double t, double T, double & g) const
{
  double A = A_->value(T);
  double n = n_->value(T);
  double beta = beta_->value(T);

  g = A * pow(sinh(beta * seq / n), n) * exp(-Q_ / (R_ * T));

  return 0;
}

int BlackburnMinimumCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
{
  double A = A_->value(T);
  double n = n_->value(T);
  double beta = beta_->value(T);

  dg = A * beta * exp(-Q_ / (R_ * T)) * cosh(beta * seq / n) * 
      pow(sinh(beta * seq / n), n - 1.0);

  return 0;
}

int BlackburnMinimumCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
{
  dg = 0.0;
  return 0;
}

int BlackburnMinimumCreep::dg_dT(double seq, double eeq, double t, double T, double & dg) const
{
  double A = A_->value(T);
  double n = n_->value(T);
  double beta = beta_->value(T);

  dg = A * pow(sinh(beta * seq / n), n) * exp(-Q_ / (R_ * T)) * Q_ / (R_ * T * T);
  return 0;
}

// Implementation of the Swindeman minimum creep rate equation
SwindemanMinimumCreep::SwindemanMinimumCreep(double C, double n, double V,
                                             double Q) :
      C_(C), n_(n), V_(V), Q_(Q)
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

  return pset;
}

std::unique_ptr<NEMLObject> SwindemanMinimumCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<SwindemanMinimumCreep>(
      params.get_parameter<double>("C"),
      params.get_parameter<double>("n"),
      params.get_parameter<double>("V"),
      params.get_parameter<double>("Q")
      ); 
}


int SwindemanMinimumCreep::g(double seq, double eeq, double t, double T, double & g) const
{
  g = C_ * pow(seq, n_) * exp(V_ * seq) * exp(-Q_/T);

  return 0;
}

int SwindemanMinimumCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
{

  dg = C_ * exp(-Q_/T) * (n_ + seq * V_) * exp(seq * V_) * pow(seq, n_ - 1.0);

  return 0;
}

int SwindemanMinimumCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
{
  dg = 0.0;
  return 0;
}

int SwindemanMinimumCreep::dg_dT(double seq, double eeq, double t, double T, double & dg) const
{
  dg = C_ * pow(seq, n_) * exp(V_ * seq) * exp(-Q_/T) * Q_ / (T * T); 
  return 0;
}

// Implementation of the mechanism switching model
RegionKMCreep::RegionKMCreep(std::vector<double> cuts, 
                             std::vector<std::shared_ptr<Interpolate>> A, 
                             std::vector<std::shared_ptr<Interpolate>> B,
                             double kboltz, double b, double eps0, 
                             std::shared_ptr<LinearElasticModel> emodel) :
    cuts_(cuts), A_(A), B_(B), kboltz_(kboltz), b_(b), eps0_(eps0), 
    b3_(pow(b,3)), emodel_(emodel)
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

  return pset;
}

std::unique_ptr<NEMLObject> RegionKMCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<RegionKMCreep>(
      params.get_parameter<std::vector<double>>("cuts"),
      params.get_object_parameter_vector<Interpolate>("A"),
      params.get_object_parameter_vector<Interpolate>("B"),
      params.get_parameter<double>("kboltz"),
      params.get_parameter<double>("b"),
      params.get_parameter<double>("eps0"),
      params.get_object_parameter<LinearElasticModel>("emodel")
      ); 
}

int RegionKMCreep::g(double seq, double eeq, double t, double T, double & g) const
{
  double A, B;
  select_region_(seq, T, A, B);
  double G = emodel_->G(T);
  double C1 = -G * b3_ / (kboltz_*T);

  g = eps0_ * exp(C1*B) * pow(seq / G, C1 * A);

  return 0;
}

int RegionKMCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
{
  double A, B;
  select_region_(seq, T, A, B);
  double G = emodel_->G(T);
  double C1 = -G * b3_ / (kboltz_*T);
  
  dg = eps0_ * exp(C1*B) * C1 * A / G * pow(seq / G, C1 * A - 1.0);

  return 0;
}

int RegionKMCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
{
  dg = 0.0;
  return 0;
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
NortonBaileyCreep::NortonBaileyCreep(std::shared_ptr<Interpolate> A,
                                     std::shared_ptr<Interpolate> m,
                                     std::shared_ptr<Interpolate> n) :
    A_(A), m_(m), n_(n)
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
  return neml::make_unique<NortonBaileyCreep>(
      params.get_object_parameter<Interpolate>("A"),
      params.get_object_parameter<Interpolate>("m"),
      params.get_object_parameter<Interpolate>("n")
      ); 
}

int NortonBaileyCreep::g(double seq, double eeq, double t, double T, double & g) const
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

  g = m * pow(A, 1.0 / m) * pow(seq, n / m) * pow(eeq, (m - 1.0) / m); 

  return 0;
}

int NortonBaileyCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
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

  dg = n * pow(A, 1.0 / m) * pow(seq, n / m - 1.0) * pow(eeq, (m - 1.0) / m);

  return 0;
}

int NortonBaileyCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
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

  dg = (m - 1) * pow(A, 1.0 / m) * pow(seq, n / m) * pow(eeq, -1.0 / m);

  return 0;
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


MukherjeeCreep::MukherjeeCreep(std::shared_ptr<LinearElasticModel> emodel,
                               double A, double n, double D0, double Q,
                               double b, double k, double R) :
    emodel_(emodel), A_(A), n_(n), D0_(D0), Q_(Q), b_(b), k_(k), R_(R)
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
  return neml::make_unique<MukherjeeCreep>(
      params.get_object_parameter<LinearElasticModel>("emodel"),
      params.get_parameter<double>("A"),
      params.get_parameter<double>("n"),
      params.get_parameter<double>("D0"),
      params.get_parameter<double>("Q"),
      params.get_parameter<double>("b"),
      params.get_parameter<double>("k"),
      params.get_parameter<double>("R")
      ); 
}

int MukherjeeCreep::g(double seq, double eeq, double t, double T, 
                      double & g) const
{
  double mu = emodel_->G(T);
  double Dv = D0_ * exp(-Q_ / (R_ * T));
  g = A_ * Dv * mu * b_ / (k_ * T) * pow(seq / mu, n_);
  return 0;
}

int MukherjeeCreep::dg_ds(double seq, double eeq, double t, double T,
                          double & dg) const
{
  double mu = emodel_->G(T);
  double Dv = D0_ * exp(-Q_ / (R_ * T));
  dg = n_ * A_ * Dv * mu * b_ / (k_ * T) * pow(seq / mu, n_ - 1.0) / mu;
  return 0;
}

int MukherjeeCreep::dg_de(double seq, double eeq, double t, double T,
                          double & dg) const
{
  dg = 0.0;
  return 0;
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

GenericCreep::GenericCreep(std::shared_ptr<Interpolate> cfn) :
    cfn_(cfn)
{

}

std::string GenericCreep::type()
{
  return "GenericCreep";
}

std::unique_ptr<NEMLObject> GenericCreep::initialize(ParameterSet & params)
{
  return neml::make_unique<GenericCreep>(
      params.get_object_parameter<Interpolate>("cfn")
      );
}

ParameterSet GenericCreep::parameters()
{
  ParameterSet pset(GenericCreep::type());

  pset.add_parameter<NEMLObject>("cfn");

  return pset;
}

int GenericCreep::g(double seq, double eeq, double t, double T, double & g) const
{
  g = exp(cfn_->value(log(seq)));
  return 0;
}

int GenericCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
{
  double f = cfn_->value(log(seq));
  double df = cfn_->derivative(log(seq));
  
  if (seq > 0.0) {
    dg = f * exp(f) * df / seq;
  }
  else {
    dg = 0.0;
  }
  return 0;
}

int GenericCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
{
  dg = 0.0;
  return 0;
}

// Setup for solve
CreepModel::CreepModel(double tol, int miter, bool verbose) :
    tol_(tol), miter_(miter), verbose_(verbose)
{

}

// Creep model default derivatives for time and temperature
int CreepModel::df_dt(const double * const s, const double * const e, double t,
                      double T, double * const df) const
{
  std::fill(df, df+6, 0.0);

  return 0;
}

int CreepModel::df_dT(const double * const s, const double * const e, double t,
                      double T, double * const df) const
{
  std::fill(df, df+6, 0.0);

  return 0;
}

// Implementation of creep model update
int CreepModel::update(const double * const s_np1, 
                       double * const e_np1, const double * const e_n,
                       double T_np1, double T_n,
                       double t_np1, double t_n,
                       double * const A_np1)
{
  // Setup the trial state
  CreepModelTrialState ts;
  int ier = make_trial_state(s_np1, e_n, T_np1, T_n, t_np1, t_n, ts);
  if (ier != SUCCESS) return ier;

  // Solve for the new creep strain
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  ier = solve(this, x, &ts, tol_, miter_, verbose_);
  if (ier != SUCCESS) return ier;
  
  // Extract
  std::copy(x, x+6, e_np1);

  // Get the tangent
  return calc_tangent_(e_np1, ts, A_np1);
}

int CreepModel::make_trial_state(const double * const s_np1, 
                                 const double * const e_n,
                                 double T_np1, double T_n,
                                 double t_np1, double t_n,
                                 CreepModelTrialState & ts) const
{
  ts.T = T_np1;
  ts.dt = t_np1 - t_n;
  ts.t = t_np1;
  std::copy(e_n, e_n+6, ts.e_n);
  std::copy(s_np1, s_np1+6, ts.s_np1);
  
  return 0;
}

// Implement the solve
size_t CreepModel::nparams() const
{
  return 6; // the creep strain
}

int CreepModel::init_x(double * const x, TrialState * ts)
{
  CreepModelTrialState * tss = static_cast<CreepModelTrialState *>(ts);

  // Just make it the previous value
  std::copy(tss->e_n, tss->e_n+6, x);
  return 0;
}

int CreepModel::RJ(const double * const x, TrialState * ts, 
                     double * const R, double * const J)
{
  CreepModelTrialState * tss = static_cast<CreepModelTrialState *>(ts);
  
  // Residual
  int ier = f(tss->s_np1, x, tss->t, tss->T, R);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<6; i++) R[i] = x[i] - tss->e_n[i] - R[i] * tss->dt;

  // Jacobian
  ier = df_de(tss->s_np1, x, tss->t, tss->T, J);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<36; i++) J[i] = -J[i] * tss->dt;
  for (int i=0; i<6; i++) J[CINDEX(i,i,6)] += 1.0;

  return 0;
}

// Helper for tangent
int CreepModel::calc_tangent_(const double * const e_np1, 
                              CreepModelTrialState & ts, double * const A_np1)
{
  int ier;
  double R[6];
  double J[36];

  ier = RJ(e_np1, &ts, R, J);
  if (ier != SUCCESS) return ier;

  ier = invert_mat(J, 6);
  if (ier != SUCCESS) return ier;

  for (int i=0; i<36; i++) J[i] = J[i] * ts.dt;

  double B[36];
  ier = df_ds(ts.s_np1, e_np1, ts.t, ts.T, B);
  if (ier != SUCCESS) return ier;

  mat_mat(6, 6, 6, J, B, A_np1);

  return 0;
}

// Implementation of J2 creep
J2CreepModel::J2CreepModel(std::shared_ptr<ScalarCreepRule> rule,
                           double tol, int miter, bool verbose) :
    CreepModel(tol, miter, verbose), rule_(rule)
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
  
  pset.add_optional_parameter<double>("tol", 1.0e-10);
  pset.add_optional_parameter<int>("miter", 25);
  pset.add_optional_parameter<bool>("verbose", false);

  return pset;
}

std::unique_ptr<NEMLObject> J2CreepModel::initialize(ParameterSet & params)
{
  return neml::make_unique<J2CreepModel>(
      params.get_object_parameter<ScalarCreepRule>("rule"),
      params.get_parameter<double>("tol"),
      params.get_parameter<int>("miter"),
      params.get_parameter<bool>("verbose")
      ); 
}

int J2CreepModel::f(const double * const s, const double * const e, double t,
                    double T, double * const f) const
{
  int ier = 0;

  // Gather the effective stresses
  double se = seq(s);
  double ee = eeq(e);

  // Get the direction
  std::copy(s, s+6, f);
  ier = sdir(f);
  if (ier != SUCCESS) return ier;

  // Get the rate
  double rate;
  ier = rule_->g(se, ee, t, T, rate);
  if (ier != SUCCESS) return ier;

  // Multiply the two together
  for (int i=0; i<6; i++) f[i] *= 3.0/2.0 * rate;

  return 0;
}

int J2CreepModel::df_ds(const double * const s, const double * const e, 
                        double t, double T, double * const df) const
{
  int ier = 0;

  // Gather the effective stresses
  double se = seq(s);
  double ee = eeq(e);
  
  // Get the direction
  double dir[6];
  std::copy(s, s+6, dir);
  ier = sdir(dir);

  // Get the rate
  double rate;
  ier = rule_->g(se, ee, t, T, rate);
  if (ier != SUCCESS) return ier;

  // Get the rate derivative
  double drate;
  ier = rule_->dg_ds(se, ee, t, T, drate);
  if (ier != SUCCESS) return ier;

  // Get our usual funny identity tensor
  double ID[36];
  std::fill(ID, ID+36, 0.0);
  for (int i=0; i<6; i++) {
    ID[CINDEX(i,i,6)] = 1.0;
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      ID[CINDEX(i,j,6)] -= 1.0 / 3.0; 
    }
  }
  // Stick on one of our 3/2 here
  for (int i=0; i<36; i++) ID[i] *= 3.0 / 2.0;

  // Begin forming outer products
  // Hack, really should figure out limit se -> 0.0
  if (se < std::numeric_limits<double>::epsilon()) {
    se = std::numeric_limits<double>::epsilon(); 
  }

  double A[36];
  ier = outer_vec(dir, 6, dir, 6, A);
  for (int i=0; i<36; i++) A[i] *= 3.0/2.0 * (drate - rate / se);
  for (int i=0; i<6; i++) A[CINDEX(i,i,6)] += rate / se;
  
  // Do the final multiplication
  mat_mat(6, 6, 6, A, ID, df);

  return 0;
}

int J2CreepModel::df_de(const double * const s, const double * const e, double t, double T, 
              double * const df) const
{
  int ier = 0;

  // Gather the effective stresses
  double se = seq(s);
  double ee = eeq(e);
  
  // Get the stress direction
  double s_dir[6];
  std::copy(s, s+6, s_dir);
  ier = sdir(s_dir);
  if (ier != SUCCESS) return ier;

  // Get the strain direction
  double e_dir[6];
  std::copy(e, e+6, e_dir);
  ier = edir(e_dir);
  if (ier != SUCCESS) return ier;

  // Get the derivative
  double drate;
  ier = rule_->dg_de(se, ee, t, T, drate);
  if (ier != SUCCESS) return ier;

  // Tack onto the direction
  for (int i=0; i<6; i++) e_dir[i] *= drate;

  // Form the final outer product
  ier = outer_vec(s_dir, 6, e_dir, 6, df);
  if (ier != SUCCESS) return ier;

  return 0;
}

int J2CreepModel::df_dt(const double * const s, const double * const e, 
                        double t, double T, double * const df) const 
{
  int ier = 0;

  // Gather the effective quantities
  double se = seq(s);
  double ee = eeq(e);
  
  // Get the stress direction
  std::copy(s, s+6, df);
  ier = sdir(df);
  if (ier != SUCCESS) return ier;

  // Get the derivative
  double drate;
  ier = rule_->dg_dt(se, ee, t, T, drate);
  if (ier != SUCCESS) return ier;

  // Multiply
  for (int i=0; i<6; i++) df[i] *= 3.0/2.0 * drate;
  
  return 0;
}

int J2CreepModel::df_dT(const double * const s, const double * const e,
                        double t, double T, double * const df) const
{
  int ier = 0;

  // Gather the effective quantities
  double se = seq(s);
  double ee = eeq(e);
  
  // Get the stress direction
  std::copy(s, s+6, df);
  ier = sdir(df);
  if (ier != SUCCESS) return ier;

  // Get the derivative
  double drate;
  ier = rule_->dg_dT(se, ee, t, T, drate);
  if (ier != SUCCESS) return ier;

  // Multiply
  for (int i=0; i<6; i++) df[i] *= 3.0/2.0 * drate;
  
  return 0;
}

// Helpers for J2 plasticity
double J2CreepModel::seq(const double * const s) const
{
  double sdev[6];
  std::copy(s, s+6, sdev);
  dev_vec(sdev);

  return sqrt(3.0 / 2.0) * norm2_vec(sdev, 6);
}

double J2CreepModel::eeq(const double * const e) const
{
  return sqrt(2.0 / 3.0) * norm2_vec(e, 6);
}

int J2CreepModel::sdir(double * const s) const
{
  int ier = 0;
  double se = seq(s);
  if (se < std::numeric_limits<double>::epsilon()) {
    std::fill(s, s+6, 0.0);
  }
  else {  
    ier = dev_vec(s);
    for (int i=0; i<6; i++) s[i] /= se;
  }
  return ier;
}

int J2CreepModel::edir(double * const e) const
{
  double ee = eeq(e);
  if (ee < std::numeric_limits<double>::epsilon()) {
    std::fill(e, e+6, 0.0);
  }
  else {
    for (int i=0; i<6; i++) e[i] /= ee;
  }
  return 0;
}

} // namespace neml
