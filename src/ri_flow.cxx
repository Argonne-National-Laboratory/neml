#include "ri_flow.h"

#include "nemlerror.h"

namespace neml {

RateIndependentFlowRule::RateIndependentFlowRule(ParameterSet & params) :
    NEMLObject(params)
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

size_t RateIndependentAssociativeFlow::nhist() const
{
  return hardening_->nhist();
}

int RateIndependentAssociativeFlow::init_hist(double * const h) const
{
  if (hardening_->nhist() != surface_->nhist()) {
    throw NEMLError("Hardening model and flow surface are not compatible");
  }

  return hardening_->init_hist(h);
}

int RateIndependentAssociativeFlow::f(const double* const s, 
                                      const double* const alpha, double T,
                                      double & fv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];

  hardening_->q(alpha, T, q);

  return surface_->f(s, q, T, fv);
}

int RateIndependentAssociativeFlow::df_ds(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  
  hardening_->q(alpha, T, q);
  
  return surface_->df_ds(s, q, T, dfv);
}

int RateIndependentAssociativeFlow::df_da(const double* const s, 
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

  return mat_vec_trans(jac, nhist(), dq, nhist(), dfv);
}

int RateIndependentAssociativeFlow::RateIndependentAssociativeFlow::g(const double * const s, 
                                      const double * const alpha, double T,
                                      double * const gv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  return surface_->df_ds(s, q, T, gv);
}

int RateIndependentAssociativeFlow::dg_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  return surface_->df_dsds(s, q, T, dgv);
}

int RateIndependentAssociativeFlow::dg_da(const double * const s, 
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

  return mat_mat(6, nhist(), nhist(), dd, jac, dgv);
}

int RateIndependentAssociativeFlow::h(const double * const s,
                                      const double * const alpha, double T,
                                      double * const hv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  return surface_->df_dq(s, q, T, hv);
}

int RateIndependentAssociativeFlow::dh_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dhv) const
{
  std::vector<double> qv(nhist());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  return surface_->df_dqds(s, q, T, dhv);
}

int RateIndependentAssociativeFlow::dh_da(const double * const s, 
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

  return mat_mat(nhist(), nhist(), nhist(), dd, jac, dhv);
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

size_t RateIndependentNonAssociativeHardening::nhist() const
{
  return hardening_->nhist();
}

int RateIndependentNonAssociativeHardening::init_hist(double * const h) const
{
  if (hardening_->ninter() != surface_->nhist()) {
    throw NEMLError("Hardening model and flow surface are not compatible");
  }

  return hardening_->init_hist(h);
}

int RateIndependentNonAssociativeHardening::f(const double* const s, 
                                      const double* const alpha, double T,
                                      double & fv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  return surface_->f(s, q, T, fv);
}

int RateIndependentNonAssociativeHardening::df_ds(const double* const s, 
                                          const double* const alpha, double T,
                                          double * const dfv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
  
  return surface_->df_ds(s, q, T, dfv);
}

int RateIndependentNonAssociativeHardening::df_da(const double* const s, 
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
  
  return mat_vec_trans(jac, nhist(), dq, hardening_->ninter(), dfv);
}

int RateIndependentNonAssociativeHardening::RateIndependentNonAssociativeHardening::g(const double * const s, 
                                      const double * const alpha, double T,
                                      double * const gv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  return surface_->df_ds(s, q, T, gv);
}

int RateIndependentNonAssociativeHardening::dg_ds(const double * const s, 
                                          const double * const alpha, double T,
                                          double * const dgv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  return surface_->df_dsds(s, q, T, dgv);
}

int RateIndependentNonAssociativeHardening::dg_da(const double * const s, 
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

  return mat_mat(6, nhist(), hardening_->ninter(), dd, jac, dgv);
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
