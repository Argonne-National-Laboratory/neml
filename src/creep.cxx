#include "creep.h"

#include "nemlmath.h"

#include <cmath>

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
PowerLawCreep::PowerLawCreep(double A, double n) :
    A_(new ConstantInterpolate(A)), n_(new ConstantInterpolate(n))
{

}

PowerLawCreep::PowerLawCreep(std::shared_ptr<Interpolate> A,
                             std::shared_ptr<Interpolate> n) :
    A_(A), n_(n)
{

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

// Implementation of Norton-Bailey creep
NortonBaileyCreep::NortonBaileyCreep(double A, double m, double n) :
    A_(new ConstantInterpolate(A)), m_(new ConstantInterpolate(m)),
    n_(new ConstantInterpolate(n))
{

}

NortonBaileyCreep::NortonBaileyCreep(std::shared_ptr<Interpolate> A,
                                     std::shared_ptr<Interpolate> m,
                                     std::shared_ptr<Interpolate> n) :
    A_(A), m_(m), n_(n)
{

}

int NortonBaileyCreep::g(double seq, double eeq, double t, double T, double & g) const
{
  double A = A_->value(T);
  double m = m_->value(T);
  double n = n_->value(T);

  g = m * pow(A, 1.0 / m) * pow(seq, n / m) * pow(eeq, (m - 1.0) / m); 

  return 0;
}

int NortonBaileyCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
{
  double A = A_->value(T);
  double m = m_->value(T);
  double n = n_->value(T);

  dg = n * pow(A, 1.0 / m) * pow(seq, n / m - 1.0) * pow(eeq, (m - 1.0) / m);

  return 0;
}

int NortonBaileyCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
{
  double A = A_->value(T);
  double m = m_->value(T);
  double n = n_->value(T);

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
                       double * const A_np1) const
{
  // Setup the trial state
  CreepModelTrialState ts;
  make_trial_state(s_np1, e_n, T_np1, T_n, t_np1, t_n, ts);

  return 0;
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
  std::copy(x, x+6, tss->e_n);
  return 0;
}

int CreepModel::RJ(const double * const x, TrialState * ts, 
                     double * const R, double * const J)
{
  CreepModelTrialState * tss = static_cast<CreepModelTrialState *>(ts);
  
  // Residual
  f(tss->s_np1, x, tss->t, tss->T, R);
  for (int i=0; i<6; i++) R[i] = x[i] - tss->e_n[i] - R[i] * tss->dt;

  // Jacobian

  return 0;
}

// Implementation of J2 creep
J2CreepModel::J2CreepModel(std::shared_ptr<ScalarCreepRule> rule) :
    rule_(rule)
{

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

  // Get the rate
  double rate;
  ier = rule_->g(se, ee, t, T, rate);

  // Multiply the two together
  for (int i=0; i<6; i++) f[i] *= 3.0/2.0 * rate;

  return ier;
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

  // Get the rate derivative
  double drate;
  ier = rule_->dg_ds(se, ee, t, T, drate);

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
  double A[36];
  ier = outer_vec(dir, 6, dir, 6, A);
  for (int i=0; i<36; i++) A[i] *= 3.0/2.0 * (drate - rate / se);
  for (int i=0; i<6; i++) A[CINDEX(i,i,6)] += rate / se;

  // Do the final multiplication
  ier = mat_mat(6, 6, 6, A, ID, df);

  return ier;
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

  // Get the strain direction
  double e_dir[6];
  std::copy(e, e+6, e_dir);
  ier = edir(e_dir);
  
  // Get the derivative
  double drate;
  ier = rule_->dg_de(se, ee, t, T, drate);

  // Tack onto the direction
  for (int i=0; i<6; i++) e_dir[i] *= drate;

  // Form the final outer product
  ier = outer_vec(s_dir, 6, e_dir, 6, df);

  return ier;
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

  // Get the derivative
  double drate;
  ier = rule_->dg_dt(se, ee, t, T, drate);

  // Multiply
  for (int i=0; i<6; i++) df[i] *= 3.0/2.0 * drate;
  
  return ier;
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

  // Get the derivative
  double drate;
  ier = rule_->dg_dT(se, ee, t, T, drate);

  // Multiply
  for (int i=0; i<6; i++) df[i] *= 3.0/2.0 * drate;
  
  return ier;
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
  double se = seq(s);
  int ier = dev_vec(s);
  for (int i=0; i<6; i++) s[i] /= se;
  return ier;
}

int J2CreepModel::edir(double * const e) const
{
  double ee = eeq(e);
  for (int i=0; i<6; i++) e[i] /= ee;
  return 0;
}

} // namespace neml
