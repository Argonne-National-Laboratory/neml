#include "visco_flow.h"

#include "nemlmath.h"

#include <cmath>
#include <iostream>

namespace neml {

// Various g(s) implementations
GPowerLaw::GPowerLaw(double n) :
    n_(n)
{

}

double GPowerLaw::g(double f) const
{
  if (f > 0.0) {
    return pow(f, n_);
  }
  else {
    return 0.0;
  }
}

double GPowerLaw::dg(double f) const
{
  if (f > 0.0) {
    return n_ * pow(f, n_ - 1.0);
  }
  else {
    return 0.0;
  }
}

double GPowerLaw::n() const
{
  return n_;
}

PerzynaFlowRule::PerzynaFlowRule(std::shared_ptr<YieldSurface> surface,
                std::shared_ptr<HardeningRule> hardening,
                std::shared_ptr<GFlow> g,
                double eta) :
    surface_(surface), hardening_(hardening), g_(g), eta_(eta)
{
  
}

size_t PerzynaFlowRule::nhist() const
{
  return hardening_->nhist();
}

int PerzynaFlowRule::init_hist(double * const h) const
{
  if (surface_->nhist() != hardening_->nhist()) {
    return INCOMPATIBLE_MODELS;
  }
  return hardening_->init_hist(h);
}

// Rate rule
int PerzynaFlowRule::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  double fv;
  surface_->f(s, q, T, fv);

  double gv = g_->g(fv);
  
  if (gv > 0.0) {
    yv = gv / eta_;
  }
  else {
    yv = 0.0;
  }

  return 0;
}

int PerzynaFlowRule::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  double fv;
  surface_->f(s, q, T, fv);

  double gv = g_->g(fv);

  std::fill(dyv, dyv + 6, 0.0);

  if (gv > 0.0) {
    double dgv = g_->dg(fv);
    surface_->df_ds(s, q, T, dyv);
    for (int i=0; i<6; i++) {
      dyv[i] *= dgv / eta_;
    }
  }
  
  return 0;
}

int PerzynaFlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  double fv;
  surface_->f(s, q, T, fv);

  double gv = g_->g(fv);
  
  std::fill(dyv, dyv + nhist(), 0.0);

  if (gv > 0.0) {
    double dgv = g_->dg(fv);

    double jac[nhist()*nhist()];
    hardening_->dq_da(alpha, T, jac);

    double rd[nhist()];
    surface_->df_dq(s, q, T, rd);

    mat_vec_trans(jac, nhist(), rd, nhist(), dyv);

    for (int i=0; i<nhist(); i++) {
      dyv[i] *= dgv / eta_;
    }
  }

  return 0;

}

// Flow rule
int PerzynaFlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->df_ds(s, q, T, gv);
}

int PerzynaFlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->df_dsds(s, q, T, dgv);
}

int PerzynaFlowRule::dg_da(const double * const s, const double * const alpha, double T,
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

// Hardening rule
int PerzynaFlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->df_dq(s, q, T, hv);
}

int PerzynaFlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  double q[nhist()];
  hardening_->q(alpha, T, q);

  return surface_->df_dqds(s, q, T, dhv);
}

int PerzynaFlowRule::dh_da(const double * const s, const double * const alpha, double T,
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

double PerzynaFlowRule::eta() const
{
  return eta_;
}

// Begin Chaboche
ConstantFluidity::ConstantFluidity(double eta) :
    eta_(eta)
{

}

double ConstantFluidity::eta(double a) const
{
  return eta_;
}

double ConstantFluidity::deta(double a) const
{
  return 0.0;
}


} // namespace neml
