#include "ri_flow.h"

#include "nemlerror.h"

namespace neml {

RateIndependentAssociativeFlow::RateIndependentAssociativeFlow(
    std::shared_ptr<YieldSurface> surface, 
    std::shared_ptr<HardeningRule> hardening) :
      surface_(surface), hardening_(hardening)
{

}

size_t RateIndependentAssociativeFlow::nhist() const
{
  return hardening_->nhist();
}

int RateIndependentAssociativeFlow::init_hist(double * const h) const
{
  if (hardening_->nhist() != surface_->nhist()) {
    return INCOMPATIBLE_MODELS;
  }

  return hardening_->init_hist(h);
}

int RateIndependentAssociativeFlow::f(const double* const s, 
                                      const double* const alpha, double T,
                                      double & fv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->f(s, q, T, fv);

}

int RateIndependentAssociativeFlow::df_ds(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);
  
  return surface_->df_ds(s, q, T, dfv);
}

int RateIndependentAssociativeFlow::df_da(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);
  
  double jac[nhist()*nhist()];
  hardening_->dq_da(alpha, T, jac);

  double dq[nhist()];
  surface_->df_dq(s, q, T, dq);

  mat_vec_trans(jac, nhist(), dq, nhist(), dfv);
}

int RateIndependentAssociativeFlow::RateIndependentAssociativeFlow::g(const double * const s, 
                                      const double * const alpha, double T,
                                      double * const gv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->df_ds(s, q, T, gv);
}

int RateIndependentAssociativeFlow::dg_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->df_dsds(s, q, T, dgv);
}

int RateIndependentAssociativeFlow::dg_da(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  double jac[nhist()*nhist()];
  hardening_->dq_da(alpha, T, jac);

  double dd[6*nhist()];
  surface_->df_dsdq(s, q, T, dd);

  mat_mat(6, nhist(), nhist(), dd, jac, dgv);

  return 0;

}

int RateIndependentAssociativeFlow::h(const double * const s,
                                      const double * const alpha, double T,
                                      double * const hv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->df_dq(s, q, T, hv);
}

int RateIndependentAssociativeFlow::dh_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->df_dqds(s, q, T, dhv);
}

int RateIndependentAssociativeFlow::dh_da(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  double jac[nhist()*nhist()];
  hardening_->dq_da(alpha, T, jac);

  double dd[nhist()*nhist()];
  surface_->df_dqdq(s, q, T, dd);

  mat_mat(nhist(), nhist(), nhist(), dd, jac, dhv);

  return 0;
}

} // namespace neml
