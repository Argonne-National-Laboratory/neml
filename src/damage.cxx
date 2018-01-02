#include "damage.h"
#include "elasticity.h"

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
  
  // Call solve
  double x[nparams()];
  int ier;
  ier = solve(this, x, &tss, tol_, miter_, verbose_);
  if (ier != SUCCESS) return ier;
  
  // Do actual stress update
  double s_np1_p[6];
  double s_n_p[6];
  double A_p[36];
  double w_n = h_n[0];
  std::copy(s_n, s_n+6, s_n_p);
  for (int i=0; i<6; i++) s_n_p[i] /= (1-w_n);
  base_->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1_p, s_n_p,
                   &h_np1[1], &h_n[1], A_p, u_np1, u_n, p_np1, p_n);
  h_np1[0] = x[6];
  for (int i=0; i<6; i++) s_np1[i] = (1-h_np1[0]) * s_np1_p[i];
  
  // Create the tangent
  ier = tangent_(e_np1, e_n, s_np1_p, s_n_p,
                 T_np1, T_n, t_np1, t_n, 
                 x[6], h_n[0], A_p, A_np1);
  if (ier != SUCCESS) return ier;

  return 0;
}

int NEMLScalarDamagedModel_sd::elastic_strains(
    const double * const s_np1,
    double T_np1,
    double * const e_np1) const
{
  // ARG this is going to require an interface change
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
}

int NEMLScalarDamagedModel_sd::RJ(const double * const x, TrialState * ts, 
                                  double * const R, double * const J)
{
  SDTrialState * tss = static_cast<SDTrialState *>(ts);
  const double * s_curr = x;
  double w_curr = x[6];
  double s_prime_curr[6];
  for (int i=0; i<6; i++)  s_prime_curr[i] = s_curr[i] / w_curr;

  double s_np1[6];
  double * h_np1 = new double [base_->nstore()];
  double A_np1[36];
  double u_np1;
  double p_np1;

  int res;

  double s_prime_n[6];
  std::copy(tss->s_n, tss->s_n+6, s_prime_n);
  for (int i=0; i<6; i++) s_prime_n[i] /= (1-tss->w_n);
  
  base_->update_sd(tss->e_np1, tss->e_n, tss->T_np1, tss->T_n,
                   tss->t_np1, tss->t_n, s_np1, s_prime_n,
                   h_np1, &tss->h_n[0],
                   A_np1, u_np1, tss->u_n, p_np1, tss->p_n);
  for (int i=0; i<6; i++) R[i] = s_curr[i] - (1-w_curr) * s_np1[i];

  double w_np1;
  damage(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, &w_np1);
  R[6] = w_curr - w_np1;

  std::fill(J, J+49, 0.0);
  for (int i=0; i<6; i++) {
    J[CINDEX(i,i,7)] = 1.0; 
  }
  for (int i=0; i<6; i++) {
    J[CINDEX(i,6,7)] = s_np1[i];
  }

  double ws[6];
  ddamage_ds(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, ws);
  for (int i=0; i<6; i++) {
    J[CINDEX(6,i,7)] = ws[i] / (1 - w_curr); 
  }
  
  double ww;
  ddamage_dd(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, &ww);
  
  J[CINDEX(6,6,7)] = 1.0 - ww;
  
  delete [] h_np1;

  return 0;
}

int NEMLScalarDamagedModel_sd::tangent_(
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n, double t_np1, double t_n,
    double w_np1, double w_n,
    const double * const A_prime, double * const A)
{
  double de[6];
  ddamage_de(w_np1, w_n, e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n, de);
  double ds[6];
  ddamage_ds(w_np1, w_n, e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n, ds);
  double dw;
  ddamage_dd(w_np1, w_n, e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n, &dw);

  double dde[6];
  mat_vec_trans(A_prime, 6, ds, 6, dde);
  for (int i=0; i<6; i++) dde[i] = (de[i] + 1.0 / (1 - w_np1) * dde[i]) / 
    (1.0 - dw);
  
  std::copy(A_prime, A_prime + 36, A);
  outer_update_minus(dde, 6, s_np1, 6, A);

  return 0;
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

int NEMLStandardScalarDamagedModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{


}

int NEMLStandardScalarDamagedModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{

}

int NEMLStandardScalarDamagedModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{

}

int NEMLStandardScalarDamagedModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{

}

}
