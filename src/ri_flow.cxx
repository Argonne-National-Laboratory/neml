#include "ri_flow.h"

#include "nemlerror.h"

namespace neml {

RateIndependentAssociativeFlow::RateIndependentAssociativeFlow(
    std::shared_ptr<YieldSurface> surface, 
    std::shared_ptr<HardeningRule> hardening) :
      surface_(surface), hardening_(hardening)
{

}

std::string RateIndependentAssociativeFlow::type()
{
  return "RateIndependentAssociativeFlow";
}

ParameterSet RateIndependentAssociativeFlow::parameters()
{
  ParameterSet pset(RateIndependentAssociativeFlow::type());

  pset.add_parameter<NEMLObject>("surface");
  pset.add_parameter<NEMLObject>("hardening");

  return pset;
}

std::shared_ptr<NEMLObject> RateIndependentAssociativeFlow::initialize(ParameterSet & params)
{
  return std::make_shared<RateIndependentAssociativeFlow>(
      params.get_object_parameter<YieldSurface>("surface"),
      params.get_object_parameter<HardeningRule>("hardening")
      ); 
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

  return 0;
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



RateIndependentNonAssociativeHardening::RateIndependentNonAssociativeHardening(
    std::shared_ptr<YieldSurface> surface, 
    std::shared_ptr<NonAssociativeHardening> hardening) :
      surface_(surface), hardening_(hardening)
{

}

std::string RateIndependentNonAssociativeHardening::type()
{
  return "RateIndependentNonAssociativeHardening";
}

ParameterSet RateIndependentNonAssociativeHardening::parameters()
{
  ParameterSet pset(RateIndependentNonAssociativeHardening::type());

  pset.add_parameter<NEMLObject>("surface");
  pset.add_parameter<NEMLObject>("hardening");

  return pset;
}

std::shared_ptr<NEMLObject> RateIndependentNonAssociativeHardening::initialize(ParameterSet & params)
{
  return std::make_shared<RateIndependentNonAssociativeHardening>(
      params.get_object_parameter<YieldSurface>("surface"),
      params.get_object_parameter<NonAssociativeHardening>("hardening")
      ); 
}

size_t RateIndependentNonAssociativeHardening::nhist() const
{
  return hardening_->nhist();
}

int RateIndependentNonAssociativeHardening::init_hist(double * const h) const
{
  if (hardening_->ninter() != surface_->nhist()) {
    return INCOMPATIBLE_MODELS;
  }

  return hardening_->init_hist(h);
}

int RateIndependentNonAssociativeHardening::f(const double* const s, 
                                      const double* const alpha, double T,
                                      double & fv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);

  return surface_->f(s, q, T, fv);

}

int RateIndependentNonAssociativeHardening::df_ds(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);
  
  return surface_->df_ds(s, q, T, dfv);
}

int RateIndependentNonAssociativeHardening::df_da(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);
  
  double jac[(hardening_->ninter())*nhist()];
  hardening_->dq_da(alpha, T, jac);

  double dq[hardening_->ninter()];
  surface_->df_dq(s, q, T, dq);
  
  // May be wrong on the ns
  mat_vec_trans(jac, nhist(), dq, hardening_->ninter(), dfv);

  return 0;
}

int RateIndependentNonAssociativeHardening::RateIndependentNonAssociativeHardening::g(const double * const s, 
                                      const double * const alpha, double T,
                                      double * const gv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);

  return surface_->df_ds(s, q, T, gv);
}

int RateIndependentNonAssociativeHardening::dg_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);

  return surface_->df_dsds(s, q, T, dgv);
}

int RateIndependentNonAssociativeHardening::dg_da(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);

  double jac[(hardening_->ninter())*nhist()];
  hardening_->dq_da(alpha, T, jac);

  double dd[6*hardening_->ninter()];
  surface_->df_dsdq(s, q, T, dd);

  mat_mat(6, nhist(), hardening_->ninter(), dd, jac, dgv);

  return 0;

}

int RateIndependentNonAssociativeHardening::h(const double * const s,
                                      const double * const alpha, double T,
                                      double * const hv) const
{
  return hardening_->h(s, alpha, T, hv);
}

int RateIndependentNonAssociativeHardening::dh_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  return hardening_->dh_ds(s, alpha, T, dhv);
}

int RateIndependentNonAssociativeHardening::dh_da(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  return hardening_->dh_da(s, alpha, T, dhv);
}

} // namespace neml
