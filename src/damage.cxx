#include "damage.h"

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
    std::shared_ptr<Interpolate> alpha) :
      NEMLDamagedModel_sd(base, alpha)
{

}

NEMLScalarDamagedModel_sd::NEMLScalarDamagedModel_sd(
    std::shared_ptr<NEMLModel_sd> base,
    double alpha) :
      NEMLDamagedModel_sd(base, alpha)
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
  for (int i=0; i<6; i++) s_prime_n[i] /= tss->w_n;
  
  base_->update_sd(tss->e_np1, tss->e_n, tss->T_np1, tss->T_n,
                   tss->t_np1, tss->t_n, s_np1, s_prime_n,
                   h_np1, &tss->h_n[0],
                   A_np1, u_np1, tss->u_n, p_np1, tss->p_n);
  for (int i=0; i<6; i++) R[i] = s_curr[i] - (1-w_curr) * s_np1[i];

  double w_np1;
  damage(w_curr, tss->e_np1, s_prime_curr, &w_np1);
  R[6] = w_curr - w_np1;

  std::fill(J, J+49, 0.0);
  for (int i=0; i<6; i++) {
    J[CINDEX(i,i,7)] = 1.0; 
  }
  for (int i=0; i<6; i++) {
    J[CINDEX(i,6,7)] = s_np1[i];
  }

  double ws[6];
  ddamage_ds(w_curr, tss->e_np1, s_prime_curr, ws);
  for (int i=0; i<6; i++) {
    J[CINDEX(6,i,7)] = ws[i] / (1 - w_curr); 
  }
  
  double ww;
  ddamage_dd(w_curr, tss->e_np1, s_prime_curr, &ww);
  
  J[CINDEX(6,6,7)] = 1.0 - ww;
  
  delete [] h_np1;

  return 0;
}

}
