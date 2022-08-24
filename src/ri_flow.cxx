#include "ri_flow.h"

#include "nemlerror.h"

namespace neml {

RateIndependentFlowRule::RateIndependentFlowRule(ParameterSet & params) :
    HistoryNEMLObject(params)
{

}

RateIndependentAssociativeFlow::RateIndependentAssociativeFlow(ParameterSet &
                                                               params) :
    RateIndependentFlowRule(params),
    surface_(params.get_object_parameter<YieldSurface>("surface")), 
    hardening_(params.get_object_parameter<HardeningRule>("hardening"))
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

std::unique_ptr<NEMLObject> RateIndependentAssociativeFlow::initialize(ParameterSet & params)
{
  return neml::make_unique<RateIndependentAssociativeFlow>(params); 
}

void RateIndependentAssociativeFlow::populate_hist(History & hist) const
{
  if (hardening_->nhist() != surface_->nhist()) {
    throw NEMLError("Hardening model and flow surface are not compatible");
  }
  hardening_->set_variable_prefix(get_variable_prefix());
  hardening_->populate_hist(hist);
}

void RateIndependentAssociativeFlow::init_hist(History & hist) const
{
  hardening_->init_hist(hist);
}

void RateIndependentAssociativeFlow::f(const double* const s, 
                                      const double* const alpha, double T,
                                      double & fv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];

  hardening_->q(alpha, T, q);

  surface_->f(s, q, T, fv);
}

void RateIndependentAssociativeFlow::df_ds(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  
  hardening_->q(alpha, T, q);
  
  surface_->df_ds(s, q, T, dfv);
}

void RateIndependentAssociativeFlow::df_da(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
  
  std::vector<double> jacv(nhist() * nhist());
  double * jac = &jacv[0];
  hardening_->dq_da(alpha, T, jac);
  
  std::vector<double> dqv(nhist());
  double * dq = &dqv[0];
  surface_->df_dq(s, q, T, dq);

  mat_vec_trans(jac, nhist(), dq, nhist(), dfv);
}

void RateIndependentAssociativeFlow::RateIndependentAssociativeFlow::g(const double * const s, 
                                      const double * const alpha, double T,
                                      double * const gv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_ds(s, q, T, gv);
}

void RateIndependentAssociativeFlow::dg_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_dsds(s, q, T, dgv);
}

void RateIndependentAssociativeFlow::dg_da(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
  
  std::vector<double> jacv(nhist() * nhist());
  double * jac = &jacv[0];
  hardening_->dq_da(alpha, T, jac);
  
  std::vector<double> ddv(6 * nhist());
  double * dd = &ddv[0];
  surface_->df_dsdq(s, q, T, dd);

  mat_mat(6, nhist(), nhist(), dd, jac, dgv);
}

void RateIndependentAssociativeFlow::h(const double * const s,
                                      const double * const alpha, double T,
                                      double * const hv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_dq(s, q, T, hv);
}

void RateIndependentAssociativeFlow::dh_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_dqds(s, q, T, dhv);
}

void RateIndependentAssociativeFlow::dh_da(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  std::vector<double> jacv(nhist() * nhist());
  double * jac = &jacv[0];
  hardening_->dq_da(alpha, T, jac);

  std::vector<double> ddv(nhist() * nhist());
  double * dd = &ddv[0];
  surface_->df_dqdq(s, q, T, dd);

  mat_mat(nhist(), nhist(), nhist(), dd, jac, dhv);
}



RateIndependentNonAssociativeHardening::RateIndependentNonAssociativeHardening(
    ParameterSet & params) :
      RateIndependentFlowRule(params),
      surface_(params.get_object_parameter<YieldSurface>("surface")),
      hardening_(params.get_object_parameter<NonAssociativeHardening>("hardening"))
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

std::unique_ptr<NEMLObject> RateIndependentNonAssociativeHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<RateIndependentNonAssociativeHardening>(params); 
}

void RateIndependentNonAssociativeHardening::populate_hist(History & hist) const
{
  if (hardening_->ninter() != surface_->nhist()) {
    throw NEMLError("Hardening model and flow surface are not compatible");
  }
  hardening_->set_variable_prefix(get_variable_prefix());
  hardening_->populate_hist(hist);
}

void RateIndependentNonAssociativeHardening::init_hist(History & hist) const
{
  hardening_->init_hist(hist);
}

void RateIndependentNonAssociativeHardening::f(const double* const s, 
                                      const double* const alpha, double T,
                                      double & fv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->f(s, q, T, fv);
}

void RateIndependentNonAssociativeHardening::df_ds(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
  
  surface_->df_ds(s, q, T, dfv);
}

void RateIndependentNonAssociativeHardening::df_da(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
 
  std::vector<double> jacv(hardening_->ninter() * nhist());
  double * jac = &jacv[0];
  hardening_->dq_da(alpha, T, jac);
  
  std::vector<double> dqv(hardening_->ninter());
  double * dq = &dqv[0];
  surface_->df_dq(s, q, T, dq);
  
  mat_vec_trans(jac, nhist(), dq, hardening_->ninter(), dfv);
}

void RateIndependentNonAssociativeHardening::RateIndependentNonAssociativeHardening::g(const double * const s, 
                                      const double * const alpha, double T,
                                      double * const gv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_ds(s, q, T, gv);
}

void RateIndependentNonAssociativeHardening::dg_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_dsds(s, q, T, dgv);
}

void RateIndependentNonAssociativeHardening::dg_da(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  std::vector<double> jacv(hardening_->ninter() * nhist());
  double * jac = &jacv[0];
  hardening_->dq_da(alpha, T, jac);
  
  std::vector<double> ddv(6 * hardening_->ninter());
  double * dd = &ddv[0];
  surface_->df_dsdq(s, q, T, dd);

  mat_mat(6, nhist(), hardening_->ninter(), dd, jac, dgv);
}

void RateIndependentNonAssociativeHardening::h(const double * const s,
                                      const double * const alpha, double T,
                                      double * const hv) const
{
  hardening_->h(s, alpha, T, hv);
}

void RateIndependentNonAssociativeHardening::dh_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  hardening_->dh_ds(s, alpha, T, dhv);
}

void RateIndependentNonAssociativeHardening::dh_da(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  hardening_->dh_da(s, alpha, T, dhv);
}

} // namespace neml
