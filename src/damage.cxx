#include "damage.h"
#include "elasticity.h"

#include <cmath>

namespace neml {

NEMLDamagedModel_sd::NEMLDamagedModel_sd(std::shared_ptr<NEMLModel_sd> base, 
                                         std::shared_ptr<Interpolate> alpha) :
    NEMLModel_sd(alpha), base_(base)
{

}

NEMLDamagedModel_sd::NEMLDamagedModel_sd(std::shared_ptr<NEMLModel_sd> base,
                                         double alpha) :
    NEMLModel_sd(alpha), base_(base)
{

}

size_t NEMLDamagedModel_sd::nhist() const
{
  return ndamage() + base_->nstore();
}

int NEMLDamagedModel_sd::init_hist(double * const hist) const
{
  int res = init_damage(hist);
  if (res != 0) return res;
  return base_->init_store(&hist[ndamage()]);
}

NEMLScalarDamagedModel_sd::NEMLScalarDamagedModel_sd(
    std::shared_ptr<NEMLModel_sd> base, 
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter, bool verbose) :
      NEMLDamagedModel_sd(base, alpha), tol_(tol), miter_(miter),
      verbose_(verbose)
{

}

NEMLScalarDamagedModel_sd::NEMLScalarDamagedModel_sd(
    std::shared_ptr<NEMLModel_sd> base,
    double alpha, double tol, int miter, bool verbose) :
      NEMLDamagedModel_sd(base, alpha), tol_(tol),
      miter_(miter), verbose_(verbose)
{

}

int NEMLScalarDamagedModel_sd::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Make trial state
  SDTrialState tss;
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n, u_n, p_n, tss);
  
  // Call solve
  double x[nparams()];
  int ier;
  ier = solve(this, x, &tss, tol_, miter_, verbose_);
  if (ier != SUCCESS) return ier;
  
  // Do actual stress update
  double s_prime_np1[6];
  double s_prime_n[6];
  double A_prime_np1[36];
  std::copy(s_n, s_n+6, s_prime_n);
  for (int i=0; i<6; i++) s_prime_n[i] /= (1-h_n[0]);

  ier = base_->update_sd(e_np1, e_n, T_np1, T_n,
                   t_np1, t_n, s_prime_np1, s_prime_n,
                   &h_np1[1], &h_n[1],
                   A_prime_np1, u_np1, u_n, p_np1, p_n);
  if (ier != SUCCESS) return ier;
  
  for (int i=0; i<6; i++) s_np1[i] = (1-x[6]) * s_prime_np1[i];
  h_np1[0] = x[6];
  
  // Create the tangent
  ier = tangent_(e_np1, e_n, s_np1, s_n,
                 T_np1, T_n, t_np1, t_n, 
                 x[6], h_n[0], A_prime_np1, A_np1);
  if (ier != SUCCESS) return ier;

  return 0;
}

size_t NEMLScalarDamagedModel_sd::ndamage() const
{
  return 1;
}

int NEMLScalarDamagedModel_sd::init_damage(double * const damage) const
{
  damage[0] = 0.0;
  return 0;
}

size_t NEMLScalarDamagedModel_sd::nparams() const
{
  return 7;
}

int NEMLScalarDamagedModel_sd::init_x(double * const x, TrialState * ts)
{
  SDTrialState * tss = static_cast<SDTrialState *>(ts);
  std::copy(tss->s_n, tss->s_n+6, x);
  x[6] = tss->w_n;

  return 0;
}

int NEMLScalarDamagedModel_sd::RJ(const double * const x, TrialState * ts, 
                                  double * const R, double * const J)
{
  SDTrialState * tss = static_cast<SDTrialState *>(ts);
  const double * s_curr = x;
  double w_curr = x[6];
  double s_prime_curr[6];
  for (int i=0; i<6; i++)  s_prime_curr[i] = s_curr[i] / (1-w_curr);

  int res;
  double s_prime_np1[6];
  double s_prime_n[6];
  double A_prime_np1[36];
  double * h_np1 = new double [base_->nstore()];
  double u_np1;
  double p_np1;
  
  std::copy(tss->s_n, tss->s_n+6, s_prime_n);
  for (int i=0; i<6; i++) s_prime_n[i] /= (1-tss->w_n);

  res = base_->update_sd(tss->e_np1, tss->e_n, tss->T_np1, tss->T_n,
                   tss->t_np1, tss->t_n, s_prime_np1, s_prime_n,
                   h_np1, &tss->h_n[0],
                   A_prime_np1, u_np1, tss->u_n, p_np1, tss->p_n);
  if (res != SUCCESS) return res;
  
  for (int i=0; i<6; i++) R[i] = s_curr[i] - (1-w_curr) * s_prime_np1[i];

  double w_np1;
  res = damage(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n, tss->t_np1, tss->t_n, &w_np1);
  if (res != SUCCESS) return res;
  R[6] = w_curr - w_np1;

  std::fill(J, J+49, 0.0);
  for (int i=0; i<6; i++) {
    J[CINDEX(i,i,7)] = 1.0; 
  }
  for (int i=0; i<6; i++) {
    J[CINDEX(i,6,7)] = s_prime_np1[i];
  }

  double ws[6];
  ddamage_ds(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, ws);
  for (int i=0; i<6; i++) {
    J[CINDEX(6,i,7)] = -ws[i] / (1 - w_curr); 
  }
  
  double ww;
  ddamage_dd(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, &ww);
  
  J[CINDEX(6,6,7)] = 1.0 - ww - dot_vec(ws, s_curr, 6) / pow(1-w_curr,2.0);

  delete [] h_np1;

  return 0;
}

int NEMLScalarDamagedModel_sd::make_trial_state(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const double * const s_n, const double * const h_n,
    double u_n, double p_n,
    SDTrialState & tss)
{
  std::copy(e_np1, e_np1+6, tss.e_np1);
  std::copy(e_n, e_n+6, tss.e_n);
  tss.T_np1 = T_np1;
  tss.T_n = T_n;
  tss.t_np1 = t_np1;
  tss.t_n = t_n;
  std::copy(s_n, s_n+6, tss.s_n);
  tss.h_n.resize(base_->nstore());
  std::copy(h_n+1, h_n+base_->nstore()+1, tss.h_n.begin());
  tss.u_n = u_n;
  tss.p_n = p_n;
  tss.w_n = h_n[0];

  return 0;
}

int NEMLScalarDamagedModel_sd::tangent_(
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n, double t_np1, double t_n,
    double w_np1, double w_n, const double * const A_prime,
    double * const A)
{
  double s_prime_np1[6];
  for (int i=0; i<6; i++) s_prime_np1[i] = s_np1[i] / (1.0 - w_np1);
  double s_prime_n[6];
  for (int i=0; i<6; i++) s_prime_n[i] = s_n[i] / (1.0 - w_n);

  double dw_ds[6];
  ddamage_ds(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, dw_ds);
  double dw_de[6];
  ddamage_de(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, dw_de);
  double dw_dw;
  ddamage_dd(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, &dw_dw);

  double k1 = 1.0 - 1.0 / (1.0 - w_np1) * dot_vec(dw_ds, s_prime_np1, 6) - dw_dw;
  double B[36];
  std::fill(B, B+36, 0);
  for (int i=0; i<6; i++) B[CINDEX(i,i,6)] = 1.0;
  for (int i=0; i<6; i++) dw_ds[i] /= (k1 * (1.0 - w_np1));
  outer_update(s_prime_np1, 6, dw_ds, 6, B);
  invert_mat(B, 6);

  double C[36];
  std::copy(A_prime, A_prime+36, C);
  for (int i=0; i<36; i++) C[i] *= (1 - w_np1);
  for (int i=0; i<6; i++) dw_de[i] /= k1;
  outer_update_minus(s_prime_np1, 6, dw_de, 6, C);
  mat_mat(6, 6, 6, B, C, A);

  return 0;
}


CombinedDamageModel_sd::CombinedDamageModel_sd(
    std::vector<std::shared_ptr<NEMLScalarDamagedModel_sd>> models,
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter, bool verbose) :
      NEMLScalarDamagedModel_sd(base, alpha, tol, miter, verbose),
      models_(models)
{


}

CombinedDamageModel_sd::CombinedDamageModel_sd(
    std::vector<std::shared_ptr<NEMLScalarDamagedModel_sd>> models,
    std::shared_ptr<NEMLModel_sd> base,
    double alpha, double tol, int miter, bool verbose) :
      NEMLScalarDamagedModel_sd(base, alpha, tol, miter, verbose),
      models_(models)
{

}

int CombinedDamageModel_sd::elastic_strains(
    const double * const s_np1,
    double T_np1, const double * const h_np1,
    double * const e_np1) const
{
  return base_->elastic_strains(s_np1, T_np1, h_np1, e_np1);
}

int CombinedDamageModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  *dd = d_n;
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    double di;
    (*it)->damage(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, &di);
    *dd += (di - d_n);
  }
  return 0;
}

int CombinedDamageModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  *dd = 0.0;
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    double di;
    (*it)->ddamage_dd(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, &di);
    *dd += di;
  }
}

int CombinedDamageModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    double di[6];
    (*it)->ddamage_de(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, di);
    for (int i=0; i<6; i++) dd[i] += di[i];
  }
}

int CombinedDamageModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    double di[6];
    (*it)->ddamage_ds(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, di);
    for (int i=0; i<6; i++) dd[i] += di[i];
  }
}

ClassicalCreepDamageModel_sd::ClassicalCreepDamageModel_sd(
    std::shared_ptr<Interpolate> A,
    std::shared_ptr<Interpolate> xi,
    std::shared_ptr<Interpolate> phi,
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter,
    bool verbose) :
      NEMLScalarDamagedModel_sd(base, alpha, tol, miter, verbose),
      A_(A), xi_(xi), phi_(phi)
{


}

ClassicalCreepDamageModel_sd::ClassicalCreepDamageModel_sd(
    double A, double xi, double phi,
    std::shared_ptr<NEMLModel_sd> base,
    double alpha,
    double tol, int miter,
    bool verbose) :
      NEMLScalarDamagedModel_sd(base, alpha, tol, miter, verbose),
      A_(new ConstantInterpolate(A)), 
      xi_(new ConstantInterpolate(xi)), 
      phi_(new ConstantInterpolate(phi)) 
{

}

int ClassicalCreepDamageModel_sd::elastic_strains(const double * const s_np1,
                    double T_np1, const double * const h_np1,
                    double * const e_np1) const
{
  return base_->elastic_strains(s_np1, T_np1, h_np1, e_np1);
}

int ClassicalCreepDamageModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se = this->se(s_np1);
  double dt = t_np1 - t_n;
  
  *dd = d_n + pow(se / A, xi) * pow(1.0 - d_np1, -phi) * dt;

  return 0;
}

int ClassicalCreepDamageModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se = this->se(s_np1);
  double dt = t_np1 - t_n;

  *dd = pow(se / A, xi) * phi * pow(1.0 - d_np1, -(phi + 1.0)) * dt;

  return 0;
}

int ClassicalCreepDamageModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);

  return 0;
}

int ClassicalCreepDamageModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se = this->se(s_np1);
  double dt = t_np1 - t_n;

  if (se == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  std::copy(s_np1, s_np1+6, dd);
  double s_m = (s_np1[0] + s_np1[1] + s_np1[2]) / 3.0;
  for (int i=0; i<3; i++) dd[i] -= s_m;

  double sf = 3.0 * xi / (2.0 * A * se) * pow(se/A, xi - 1.0) * 
      pow(1 - d_np1, -phi) * dt;
  for (int i=0; i<6; i++) dd[i] *= sf;
  
  return 0;
}

double ClassicalCreepDamageModel_sd::se(const double * const s) const
{
  return sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) + 
               pow(s[2] - s[0], 2.0) + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + 
                                              pow(s[5], 2.0))) / 2.0);
}

MarkFatigueDamageModel_sd::MarkFatigueDamageModel_sd(
    std::shared_ptr<Interpolate> C,
    std::shared_ptr<Interpolate> m,
    std::shared_ptr<Interpolate> n,
    std::shared_ptr<Interpolate> falpha,
    std::shared_ptr<Interpolate> fbeta,
    std::shared_ptr<Interpolate> rate0,
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter,
    bool verbose) :
      NEMLScalarDamagedModel_sd(base, alpha, tol, miter, verbose),
      C_(C), m_(m), n_(n), falpha_(falpha), fbeta_(fbeta), rate0_(rate0)
{


}

MarkFatigueDamageModel_sd::MarkFatigueDamageModel_sd(
    double C, double m, double n,
    double falpha, double fbeta, double rate0,
    std::shared_ptr<NEMLModel_sd> base,
    double alpha,
    double tol, int miter,
    bool verbose) :
      NEMLScalarDamagedModel_sd(base, alpha, tol, miter, verbose),
      C_(new ConstantInterpolate(C)), 
      m_(new ConstantInterpolate(m)), 
      n_(new ConstantInterpolate(n)), 
      falpha_(new ConstantInterpolate(falpha)), 
      fbeta_(new ConstantInterpolate(fbeta)),
      rate0_(new ConstantInterpolate(rate0))
{

}

int MarkFatigueDamageModel_sd::elastic_strains(const double * const s_np1,
                    double T_np1, const double * const h_np1,
                    double * const e_np1) const
{
  return base_->elastic_strains(s_np1, T_np1, h_np1, e_np1);
}

int MarkFatigueDamageModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double C = C_->value(T_np1);
  double m = m_->value(T_np1);
  double n = n_->value(T_np1);
  double fa = falpha_->value(T_np1);
  double fb = fbeta_->value(T_np1);
  double r0 = rate0_->value(T_np1);

  double de[6];
  for (int i=0; i<6; i++) de[i] = e_np1[i] - e_n[i];

  double seq = se_(s_np1);
  double eeq = ee_(de);

  if ((seq == 0) || (eeq == 0)) {
    *dd = d_n;
    return 0;
  }

  double dt = t_np1 - t_n;
  if (dt == 0.0) {
    *dd = d_n;
    return 0;
  }

  *dd = d_n + C * beta_fn_(d_np1, fa, fb, r0) * pow(seq, n) * pow(eeq/dt, m) * dt;

  return 0;
}

int MarkFatigueDamageModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double C = C_->value(T_np1);
  double m = m_->value(T_np1);
  double n = n_->value(T_np1);
  double fa = falpha_->value(T_np1);
  double fb = fbeta_->value(T_np1);
  double r0 = rate0_->value(T_np1);

  double de[6];
  for (int i=0; i<6; i++) de[i] = e_np1[i] - e_n[i];

  double seq = se_(s_np1);
  double eeq = ee_(de);

  if ((seq == 0) || (eeq == 0)) {
    *dd = 0.0;
    return 0;
  }

  double dt = t_np1 - t_n;
  if (dt == 0.0) {
    *dd = 0.0;
    return 0;
  }

  *dd = C * d_beta_fn_(d_np1, fa, fb, r0) * pow(seq, n) * pow(eeq/dt, m) * dt;

  return 0;
}

int MarkFatigueDamageModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double C = C_->value(T_np1);
  double m = m_->value(T_np1);
  double n = n_->value(T_np1);
  double fa = falpha_->value(T_np1);
  double fb = fbeta_->value(T_np1);
  double r0 = rate0_->value(T_np1);

  double de[6];
  for (int i=0; i<6; i++) de[i] = e_np1[i] - e_n[i];

  double seq = se_(s_np1);
  double eeq = ee_(de);

  if ((seq == 0) || (eeq == 0)) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  double dt = t_np1 - t_n;
  if (dt == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  dee_(de, dd);
  double fact = C * beta_fn_(d_np1, fa, fb, r0) * pow(seq, n) * m * pow(eeq/dt, m-1);

  for (int i=0; i<6; i++) dd[i] *= fact;

  return 0;
}

int MarkFatigueDamageModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double C = C_->value(T_np1);
  double m = m_->value(T_np1);
  double n = n_->value(T_np1);
  double fa = falpha_->value(T_np1);
  double fb = fbeta_->value(T_np1);
  double r0 = rate0_->value(T_np1);

  double de[6];
  for (int i=0; i<6; i++) de[i] = e_np1[i] - e_n[i];

  double seq = se_(s_np1);
  double eeq = ee_(de);

  if ((seq == 0) || (eeq == 0)) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  double dt = t_np1 - t_n;
  if (dt == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  dse_(s_np1, dd);
  double fact = C * beta_fn_(d_np1, fa, fb, r0) * n * pow(seq, n-1) * pow(eeq/dt, m) * dt;

  for (int i=0; i<6; i++) dd[i] *= fact;

  return 0;

}

double MarkFatigueDamageModel_sd::beta_fn_(double x, double a, double b, double r0) const
{
  double nz = tgamma(a+b) / (tgamma(a) + tgamma(b));
  return nz * pow(x, a-1) * pow(1-x, b-1) + r0 * (1-x);
}

double MarkFatigueDamageModel_sd::d_beta_fn_(double x, double a, double b, double r0) const
{
  double nz = tgamma(a+b) / (tgamma(a) + tgamma(b));
  return nz * (
      (a-1) * pow(x, a-2) * pow(1-x, b-1) - 
      (b-1) * pow(x, a-1) * pow(1-x, b-2)) - r0;
}

double MarkFatigueDamageModel_sd::se_(const double * const s) const
{
  double sdev[6];
  std::copy(s, s+6, sdev);
  dev_vec(sdev);
  return sqrt(3.0/2.0 * dot_vec(sdev, sdev, 6));
}

void MarkFatigueDamageModel_sd::dse_(const double * const s, double * const ds) const
{
  std::copy(s, s+6, ds);
  dev_vec(ds);
  double f = dot_vec(ds, ds, 6);
  if (f == 0.0) {
    std::fill(ds, ds+6, 0.0);
    return;
  }
  else {
    double fact = sqrt(3.0/2.0) / sqrt(f);
    for (int i = 0; i < 6; i++) ds[i] *= fact;
    return;
  }
}

double MarkFatigueDamageModel_sd::ee_(const double * const e) const
{
  double edev[6];
  std::copy(e, e+6, edev);
  dev_vec(edev);
  return sqrt(2.0/3.0 * dot_vec(edev, edev, 6));
}

void MarkFatigueDamageModel_sd::dee_(const double * const e, double * const de) const
{
  std::copy(e, e+6, de);
  dev_vec(de);
  double f = dot_vec(de, de, 6);
  if (f == 0.0) {
    std::fill(de, de+6, 0.0);
    return;
  }
  else {
    double fact = sqrt(2.0/3.0) / sqrt(f);
    for (int i = 0; i < 6; i++) de[i] *= fact;
    return;
  }
}

NEMLStandardScalarDamagedModel_sd::NEMLStandardScalarDamagedModel_sd(
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<LinearElasticModel> emodel,
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter, bool verbose) :
      NEMLScalarDamagedModel_sd(base, alpha, tol, miter, verbose), 
      emodel_(emodel)
{

}

NEMLStandardScalarDamagedModel_sd::NEMLStandardScalarDamagedModel_sd(
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<LinearElasticModel> emodel,
    double alpha, double tol, int miter, bool verbose) :
      NEMLScalarDamagedModel_sd(base, alpha, tol, miter, verbose),
      emodel_(emodel)
{

}

int NEMLStandardScalarDamagedModel_sd::elastic_strains(
    const double * const s_np1, double T_np1,
    const double * const h_np1, double * const e_np1) const
{
  double Sv[36];
  emodel_->S(T_np1, Sv);
  for (int i=0; i<36; i++) Sv[i] /= (1 - h_np1[0]);
  return mat_vec(Sv, 6, s_np1, 6, e_np1);
}

int NEMLStandardScalarDamagedModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  int ier = f(s_np1, d_np1, T_np1, fval);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);
  *dd = d_n + fval * deps;

  return 0;
}

int NEMLStandardScalarDamagedModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double df;
  int ier = df_dd(s_np1, d_np1, T_np1, df);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  *dd = df * deps;

  return 0;
}

int NEMLStandardScalarDamagedModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  int ier = f(s_np1, d_np1, T_np1, fval);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  if (deps == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }
  
  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] - s_n[i];
    de[i] = e_np1[i] - e_n[i];
  }

  double S[36];
  emodel_->S(T_np1, S);

  double dee[36];
  mat_vec(S, 6, ds, 6, dee);

  for (int i=0; i<6; i++) {
    dd[i] = (2.0 * fval) / (3.0 * deps) * (de[i] - dee[i]); 
  }

  return 0;
}

int NEMLStandardScalarDamagedModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  int ier = f(s_np1, d_np1, T_np1, fval);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  if (deps == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] - s_n[i];
    de[i] = e_np1[i] - e_n[i];
  }

  double S[36];
  emodel_->S(T_np1, S);

  double dee[36];
  mat_vec(S, 6, ds, 6, dee);

  double v1[6];
  for (int i=0; i<6; i++) {
    v1[i] = (2.0 * fval) / (3.0 * deps) * (dee[i] - de[i]); 
  }

  mat_vec(S, 6, v1, 6, dd);

  double dds[6];
  df_ds(s_np1, d_np1, T_np1, dds);

  for (int i=0; i<6; i++) {
    dd[i] = dds[i] * deps + dd[i];
  }

  return 0;
}

double NEMLStandardScalarDamagedModel_sd::dep(
    const double * const s_np1, const double * const s_n,
    const double * const e_np1, const double * const e_n,
    double T_np1) const
{
  double S[36];
  emodel_->S(T_np1, S);
  
  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] - s_n[i];
    de[i] = e_np1[i] - e_n[i];
  }

  double dee[36];
  mat_vec(S, 6, ds, 6, dee);

  double val_sq =  2.0/3.0 * (dot_vec(de, de, 6) + dot_vec(dee, dee, 6) -
                         2.0 * dot_vec(de, dee, 6));
  
  if (val_sq < 0.0) {
    return 0.0;
  }
  else {
    return sqrt(val_sq);
  }
}

NEMLPowerLawDamagedModel_sd::NEMLPowerLawDamagedModel_sd(
    std::shared_ptr<Interpolate> A, std::shared_ptr<Interpolate> a, 
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<LinearElasticModel> emodel,
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter,
    bool verbose) :
      NEMLStandardScalarDamagedModel_sd(base, emodel, alpha, tol, miter, 
                                        verbose), 
      A_(A), a_(a)
{

}

NEMLPowerLawDamagedModel_sd::NEMLPowerLawDamagedModel_sd(
    double A, double a,
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<LinearElasticModel> emodel,
    double alpha,
    double tol, int miter,
    bool verbose) :
      NEMLStandardScalarDamagedModel_sd(base, emodel, alpha, tol, miter,
                                        verbose),
      A_(new ConstantInterpolate(A)), a_(new ConstantInterpolate(a))
{

}

int NEMLPowerLawDamagedModel_sd::f(const double * const s_np1, double d_np1,
                                double T_np1, double & f) const
{
  double sev = se(s_np1);
  double A = A_->value(T_np1);
  double a = a_->value(T_np1);

  f = A * pow(sev, a);

  return 0;
}

int NEMLPowerLawDamagedModel_sd::df_ds(const double * const s_np1, double d_np1, double T_np1,
                                 double * const df) const
{
  double sev = se(s_np1);
  double A = A_->value(T_np1);
  double a = a_->value(T_np1);

  if (sev == 0.0) {
    std::fill(df, df+6, 0.0);
    return 0.0;
  }

  std::copy(s_np1, s_np1+6, df);
  double sm = (s_np1[0] + s_np1[1] + s_np1[2]) / 3.0;
  for (int i=0; i<3; i++) df[i] -= sm;
  for (int i=0; i<6; i++) df[i] *= (3.0 * A * a / 2.0 * pow(sev, a - 2.0));

  return 0;
}

int NEMLPowerLawDamagedModel_sd::df_dd(const double * const s_np1, double d_np1, double T_np1,
                                 double & df) const
{
  df = 0.0;

  return 0;
}

double NEMLPowerLawDamagedModel_sd::se(const double * const s) const
{
  return sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) + 
               pow(s[2] - s[0], 2.0) + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + 
                                              pow(s[5], 2.0))) / 2.0);
}

NEMLExponentialWorkDamagedModel_sd::NEMLExponentialWorkDamagedModel_sd(
    std::shared_ptr<Interpolate> W0, std::shared_ptr<Interpolate> k0,
    std::shared_ptr<Interpolate> af,
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<LinearElasticModel> emodel,
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter,
    bool verbose) :
      NEMLStandardScalarDamagedModel_sd(base, emodel, alpha, tol, miter, 
                                        verbose), 
      W0_(W0), k0_(k0), af_(af)
{

}

NEMLExponentialWorkDamagedModel_sd::NEMLExponentialWorkDamagedModel_sd(
    double W0, double k0, double af,
    std::shared_ptr<NEMLModel_sd> base,
    std::shared_ptr<LinearElasticModel> emodel,
    double alpha,
    double tol, int miter,
    bool verbose) :
      NEMLStandardScalarDamagedModel_sd(base, emodel, alpha, tol, miter,
                                        verbose),
      W0_(new ConstantInterpolate(W0)), k0_(new ConstantInterpolate(k0)),
      af_(new ConstantInterpolate(af))
{

}

int NEMLExponentialWorkDamagedModel_sd::f(const double * const s_np1, double d_np1,
                                double T_np1, double & f) const
{
  double sev = se(s_np1);
  double W0 = W0_->value(T_np1);
  double k0 = k0_->value(T_np1);
  double af = af_->value(T_np1);
  
  // Sign, odd case during iteration
  if ((d_np1 + k0) < 0.0) {
    f = 0.0;
  }
  else {
    f = pow(d_np1 + k0, af) / W0 * sev;
  }

  return 0;
}

int NEMLExponentialWorkDamagedModel_sd::df_ds(const double * const s_np1, double d_np1, double T_np1,
                                 double * const df) const
{
  double sev = se(s_np1);
  double W0 = W0_->value(T_np1);
  double k0 = k0_->value(T_np1);
  double af = af_->value(T_np1);

  if (sev == 0.0) {
    std::fill(df, df+6, 0.0);
    return 0;
  }

  if ((d_np1 + k0) < 0.0) {
    std::fill(df, df+6, 0.0);
    return 0;
  }

  std::copy(s_np1, s_np1+6, df);
  double sm = (s_np1[0] + s_np1[1] + s_np1[2]) / 3.0;
  for (int i=0; i<3; i++) df[i] -= sm;
  for (int i=0; i<6; i++) df[i] *= (3.0 * pow(d_np1 + k0, af) / (2.0 * sev * W0) );

  return 0;
}

int NEMLExponentialWorkDamagedModel_sd::df_dd(const double * const s_np1, double d_np1, double T_np1,
                                 double & df) const
{
  double sev = se(s_np1);
  double W0 = W0_->value(T_np1);
  double k0 = k0_->value(T_np1);
  double af = af_->value(T_np1);

  if ((d_np1 + k0) < 0.0) {
    df = 0.0;
    return 0;
  }

  df = af * pow(d_np1 + k0, af - 1.0) * sev / W0;

  return 0;
}

double NEMLExponentialWorkDamagedModel_sd::se(const double * const s) const
{
  return sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) + 
               pow(s[2] - s[0], 2.0) + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + 
                                              pow(s[5], 2.0))) / 2.0);
}

}
