#include "visco_flow.h"

#include "nemlmath.h"

#include <cmath>
#include <iostream>

namespace neml {

// Default implementation of flow rule wrt time
int ViscoPlasticFlowRule::g_time(const double * const s, 
                                 const double * const alpha, double T, 
                                 double * const gv) const
{
  std::fill(gv, gv+6, 0.0);
  return 0;
}

int ViscoPlasticFlowRule::dg_ds_time(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dgv) const
{
  std::fill(dgv, dgv+36, 0.0);
  return 0;
}

int ViscoPlasticFlowRule::dg_da_time(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dgv) const
{
  std::fill(dgv, dgv+6*nhist(), 0.0);
  return 0;
}

// Default implementation of flow rule wrt temperature
int ViscoPlasticFlowRule::g_temp(const double * const s, 
                                 const double * const alpha, double T, 
                                 double * const gv) const
{
  std::fill(gv, gv+6, 0.0);
  return 0;
}

int ViscoPlasticFlowRule::dg_ds_temp(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dgv) const
{
  std::fill(dgv, dgv+36, 0.0);
  return 0;
}

int ViscoPlasticFlowRule::dg_da_temp(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dgv) const
{
  std::fill(dgv, dgv+6*nhist(), 0.0);
  return 0;
}

// Default implementation of hardening rule wrt time
int ViscoPlasticFlowRule::h_time(const double * const s, 
                                 const double * const alpha, double T, 
                                 double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
  return 0;
}

int ViscoPlasticFlowRule::dh_ds_time(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dhv) const
{
  std::fill(dhv, dhv+6*nhist(), 0.0);
  return 0;
}

int ViscoPlasticFlowRule::dh_da_time(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);
  return 0;
}

// Default implementation of hardening rule wrt temperature
int ViscoPlasticFlowRule::h_temp(const double * const s, 
                                 const double * const alpha, double T, 
                                 double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
  return 0;
}

int ViscoPlasticFlowRule::dh_ds_temp(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dhv) const
{
  std::fill(dhv, dhv+6*nhist(), 0.0);
  return 0;
}

int ViscoPlasticFlowRule::dh_da_temp(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);
  return 0;
}

// Various g(s) implementations
GPowerLaw::GPowerLaw(double n) :
    n_(new ConstantInterpolate(n))
{

}

GPowerLaw::GPowerLaw(std::shared_ptr<Interpolate> n) :
    n_(n)
{

}

std::string GPowerLaw::type()
{
  return "GPowerLaw";
}

ParameterSet GPowerLaw::parameters()
{
  ParameterSet pset(GPowerLaw::type());

  pset.add_parameter<std::shared_ptr<Interpolate>>("n");

  return pset;
}

std::shared_ptr<NEMLObject> GPowerLaw::initialize(ParameterSet & params)
{
  return std::make_shared<GPowerLaw>(
      params.get_parameter<std::shared_ptr<Interpolate>>("n")
      ); 
}

double GPowerLaw::g(double f, double T) const
{
  if (f > 0.0) {
    return pow(f, n_->value(T));
  }
  else {
    return 0.0;
  }
}

double GPowerLaw::dg(double f, double T) const
{
  if (f > 0.0) {
    return n_->value(T) * pow(f, n_->value(T) - 1.0);
  }
  else {
    return 0.0;
  }
}

double GPowerLaw::n(double T) const
{
  return n_->value(T);
}

PerzynaFlowRule::PerzynaFlowRule(std::shared_ptr<YieldSurface> surface,
                std::shared_ptr<HardeningRule> hardening,
                std::shared_ptr<GFlow> g,
                double eta) :
    surface_(surface), hardening_(hardening), g_(g), 
    eta_(new ConstantInterpolate(eta))
{
  
}

PerzynaFlowRule::PerzynaFlowRule(std::shared_ptr<YieldSurface> surface,
                std::shared_ptr<HardeningRule> hardening,
                std::shared_ptr<GFlow> g,
                std::shared_ptr<Interpolate> eta) :
    surface_(surface), hardening_(hardening), g_(g), eta_(eta)
{
  
}

std::string PerzynaFlowRule::type()
{
  return "PerzynaFlowRule";
}

ParameterSet PerzynaFlowRule::parameters()
{
  ParameterSet pset(PerzynaFlowRule::type());

  pset.add_parameter<std::shared_ptr<NEMLObject>>("surface");
  pset.add_parameter<std::shared_ptr<NEMLObject>>("hardening");
  pset.add_parameter<std::shared_ptr<NEMLObject>>("g");
  pset.add_parameter<std::shared_ptr<Interpolate>>("eta");

  return pset;
}

std::shared_ptr<NEMLObject> PerzynaFlowRule::initialize(ParameterSet & params)
{
  return std::make_shared<PerzynaFlowRule>(
      params.get_object_parameter<YieldSurface>("surface"),
      params.get_object_parameter<HardeningRule>("hardening"),
      params.get_object_parameter<GFlow>("g"),
      params.get_parameter<std::shared_ptr<Interpolate>>("eta")
      ); 
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

  double gv = g_->g(fv, T);
  
  if (gv > 0.0) {
    yv = gv / eta_->value(T);
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

  double gv = g_->g(fv, T);

  std::fill(dyv, dyv + 6, 0.0);

  if (gv > 0.0) {
    double dgv = g_->dg(fv, T);
    surface_->df_ds(s, q, T, dyv);
    for (int i=0; i<6; i++) {
      dyv[i] *= dgv / eta_->value(T);
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

  double gv = g_->g(fv, T);
  
  std::fill(dyv, dyv + nhist(), 0.0);

  if (gv > 0.0) {
    double dgv = g_->dg(fv, T);

    double jac[nhist()*nhist()];
    hardening_->dq_da(alpha, T, jac);

    double rd[nhist()];
    surface_->df_dq(s, q, T, rd);

    mat_vec_trans(jac, nhist(), rd, nhist(), dyv);

    for (int i=0; i<nhist(); i++) {
      dyv[i] *= dgv / eta_->value(T);
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

double PerzynaFlowRule::eta(double T) const
{
  return eta_->value(T);
}

// Begin Chaboche
ConstantFluidity::ConstantFluidity(double eta) :
    eta_(new ConstantInterpolate(eta))
{

}

ConstantFluidity::ConstantFluidity(std::shared_ptr<Interpolate> eta) :
    eta_(eta)
{

}

std::string ConstantFluidity::type()
{
  return "ConstantFluidity";
}

ParameterSet ConstantFluidity::parameters()
{
  ParameterSet pset(ConstantFluidity::type());

  pset.add_parameter<std::shared_ptr<Interpolate>>("eta");

  return pset;
}

std::shared_ptr<NEMLObject> ConstantFluidity::initialize(ParameterSet & params)
{
  return std::make_shared<ConstantFluidity>(
      params.get_parameter<std::shared_ptr<Interpolate>>("eta")
      ); 
}


double ConstantFluidity::eta(double a, double T) const
{
  return eta_->value(T);
}

double ConstantFluidity::deta(double a, double T) const
{
  return 0.0;
}

SaturatingFluidity::SaturatingFluidity(double K0, double A, double b)
  : K0_(new ConstantInterpolate(K0)),
    A_(new ConstantInterpolate(A)),
    b_(new ConstantInterpolate(b))
{

}

SaturatingFluidity::SaturatingFluidity(std::shared_ptr<Interpolate> K0,
                   std::shared_ptr<Interpolate> A,
                   std::shared_ptr<Interpolate> b)
  : K0_(K0), A_(A), b_(b)
{

}

std::string SaturatingFluidity::type()
{
  return "SaturatingFluidity";
}

ParameterSet SaturatingFluidity::parameters()
{
  ParameterSet pset(SaturatingFluidity::type());

  pset.add_parameter<std::shared_ptr<Interpolate>>("K0");
  pset.add_parameter<std::shared_ptr<Interpolate>>("A");
  pset.add_parameter<std::shared_ptr<Interpolate>>("b");

  return pset;
}

std::shared_ptr<NEMLObject> SaturatingFluidity::initialize(ParameterSet & params)
{
  return std::make_shared<SaturatingFluidity>(
      params.get_parameter<std::shared_ptr<Interpolate>>("K0"),
      params.get_parameter<std::shared_ptr<Interpolate>>("A"),
      params.get_parameter<std::shared_ptr<Interpolate>>("b")
      ); 
}


double SaturatingFluidity::eta(double a, double T) const
{
  double K0 = K0_->value(T);
  double A = A_->value(T);
  double b = b_->value(T);

  return K0 + A * (1.0 - exp(-b * a));
}

double SaturatingFluidity::deta(double a, double T) const
{
  double A = A_->value(T);
  double b = b_->value(T);

  return A * b * exp(-b * a);
}

ChabocheFlowRule::ChabocheFlowRule(std::shared_ptr<YieldSurface> surface,
                                   std::shared_ptr<NonAssociativeHardening> hardening,
                                   std::shared_ptr<FluidityModel> fluidity,
                                   double n) :
    surface_(surface), hardening_(hardening), fluidity_(fluidity), 
    n_(new ConstantInterpolate(n)), recovery_(false)
{
  
}

ChabocheFlowRule::ChabocheFlowRule(std::shared_ptr<YieldSurface> surface,
                                   std::shared_ptr<NonAssociativeHardening> hardening,
                                   std::shared_ptr<FluidityModel> fluidity,
                                   std::shared_ptr<Interpolate> n) :
    surface_(surface), hardening_(hardening), fluidity_(fluidity), n_(n),
    recovery_(false)
{
  
}

std::string ChabocheFlowRule::type()
{
  return "ChabocheFlowRule";
}

ParameterSet ChabocheFlowRule::parameters()
{
  ParameterSet pset(ChabocheFlowRule::type());

  pset.add_parameter<std::shared_ptr<NEMLObject>>("surface");
  pset.add_parameter<std::shared_ptr<NEMLObject>>("hardening");
  pset.add_parameter<std::shared_ptr<NEMLObject>>("fluidity");
  pset.add_parameter<std::shared_ptr<Interpolate>>("n");

  return pset;
}

std::shared_ptr<NEMLObject> ChabocheFlowRule::initialize(ParameterSet & params)
{
  return std::make_shared<ChabocheFlowRule>(
      params.get_object_parameter<YieldSurface>("surface"),
      params.get_object_parameter<NonAssociativeHardening>("hardening"),
      params.get_object_parameter<FluidityModel>("fluidity"),
      params.get_parameter<std::shared_ptr<Interpolate>>("n")
      ); 
}

size_t ChabocheFlowRule::nhist() const
{
  return hardening_->nhist();
}

int ChabocheFlowRule::init_hist(double * const h) const
{
  if (surface_->nhist() != hardening_->ninter()) {
    return INCOMPATIBLE_MODELS;
  }
  return hardening_->init_hist(h);
}

// Rate rule
int ChabocheFlowRule::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);
  
  double fv;
  surface_->f(s, q, T, fv);

  if (fv > 0.0) {
    double eta = sqrt(2.0/3.0) * fluidity_->eta(alpha[0], T);
    yv = sqrt(3.0/2.0) * pow(fv/eta, n_->value(T));
  }
  else {
    yv = 0.0;
  }

  return 0;
}

int ChabocheFlowRule::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);
 
  double fv;
  surface_->f(s, q, T, fv);
  
  std::fill(dyv, dyv+6, 0.0);

  if (fv > 0.0) {
    surface_->df_ds(s, q, T, dyv);
    double eta = sqrt(2.0/3.0) * fluidity_->eta(alpha[0], T);
    double mv = sqrt(3.0/2.0) * pow(fv/eta, n_->value(T) - 1.0) * n_->value(T) / eta;
    for (int i=0; i<6; i++) dyv[i] *= mv;
  }

  return 0;
}

int ChabocheFlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);

  double fv;
  surface_->f(s, q, T, fv);

  std::fill(dyv, dyv + nhist(), 0.0);

  if (fv > 0.0) {
    double jac[(hardening_->ninter())*nhist()];
    hardening_->dq_da(alpha, T, jac);

    double dq[hardening_->ninter()];
    surface_->df_dq(s, q, T, dq);

    mat_vec_trans(jac, nhist(), dq, hardening_->ninter(), dyv);

    double eta = sqrt(2.0/3.0) * fluidity_->eta(alpha[0], T);
    double mv = sqrt(3.0/2.0) * pow(fv/eta, n_->value(T) - 1.0) * n_->value(T) / eta;
    for (int i=0; i<nhist(); i++) dyv[i] *= mv;

    double mv2 = -sqrt(3.0/2.0) * fv * pow(fv/eta, n_->value(T) - 1.0) * n_->value(T) / (eta * eta);
    double deta = sqrt(2.0/3.0) * fluidity_->deta(alpha[0], T);
    dyv[0] += deta * mv2;
  }

  return 0;

}

// Flow rule
int ChabocheFlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);

  return surface_->df_ds(s, q, T, gv);
}

int ChabocheFlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  double q[hardening_->ninter()];
  hardening_->q(alpha, T, q);

  return surface_->df_dsds(s, q, T, dgv);
}

int ChabocheFlowRule::dg_da(const double * const s, const double * const alpha, double T,
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

// Hardening rule
int ChabocheFlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  return hardening_->h(s, alpha, T, hv);
}

int ChabocheFlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return hardening_->dh_ds(s, alpha, T, dhv);
}

int ChabocheFlowRule::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return hardening_->dh_da(s, alpha, T, dhv);
}

// Hardening rule wrt time
int ChabocheFlowRule::h_time(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  return hardening_->h_time(s, alpha, T, hv);
}

int ChabocheFlowRule::dh_ds_time(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return hardening_->dh_ds_time(s, alpha, T, dhv);
}

int ChabocheFlowRule::dh_da_time(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return hardening_->dh_da_time(s, alpha, T, dhv);
}

// Hardening rule wrt temperature
int ChabocheFlowRule::h_temp(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  return hardening_->h_temp(s, alpha, T, hv);
}

int ChabocheFlowRule::dh_ds_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return hardening_->dh_ds_temp(s, alpha, T, dhv);
}

int ChabocheFlowRule::dh_da_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return hardening_->dh_da_temp(s, alpha, T, dhv);
}

YaguchiGr91FlowRule::YaguchiGr91FlowRule()
{

}

std::string YaguchiGr91FlowRule::type()
{
  return "YaguchiGr91FlowRule";
}

ParameterSet YaguchiGr91FlowRule::parameters()
{
  ParameterSet pset(YaguchiGr91FlowRule::type());

  return pset;
}

std::shared_ptr<NEMLObject> YaguchiGr91FlowRule::initialize(ParameterSet & params)
{
  return std::make_shared<YaguchiGr91FlowRule>(); 
}


size_t YaguchiGr91FlowRule::nhist() const
{
  // Order:
  // 0-5  X1
  // 6-11 X2
  // 12   Q
  // 13   sa
  return 14;
}

int YaguchiGr91FlowRule::init_hist(double * const h) const
{
  // This also hardcoded from the paper
  
  // Xs
  std::fill(h, h+12, 0.0);

  // Q
  h[12] = 0.0;

  // sa
  h[13] = 0.0;

  return 0;

}

// Rate rule
int YaguchiGr91FlowRule::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  double nT = n(T);
  double DT = D(T);
  double sa = alpha[13];

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  yv = (J2_(dS) - sa) / DT;
  if (yv > 0.0) {
    yv = pow(fabs(yv), nT);
  }
  else {
    yv = 0.0;
  }

  return 0;
}

int YaguchiGr91FlowRule::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  std::fill(dyv, dyv+6, 0.0);

  double yi;
  y(s, alpha, T, yi);
  double nT = n(T);
  double DT = D(T);
  double sa = alpha[13];

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);
  
  if (yi > 0.0) {
    double j2 = J2_(dS);
    double sp = (j2 - sa) / DT;
    sp = pow(fabs(sp), nT - 1.0) * nT * copysign(1.0, sp) / DT;
    dev_vec_deriv_(dS, dyv);
    for (int i=0; i<6; i++) {
      dyv[i] *= 3.0/2.0 / j2 * sp;
    }
  }
  else {
    std::fill(dyv, dyv+6, 0.0);
  }

  return 0;
}

int YaguchiGr91FlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  std::fill(dyv, dyv+nhist(), 0.0);

  // General
  double yi;
  y(s, alpha, T, yi);
  double nT = n(T);
  double DT = D(T);
  double sa = alpha[13];

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);
  double j2 = J2_(dS);
  double q = (j2 - sa) / DT;

  if (yi > 0.0) {
    // Xs
    double sp = (j2 - sa) / DT;
    sp = pow(fabs(sp), nT - 1.0) * nT * copysign(1.0, sp) / DT;

    dev_vec_deriv_(dS, &dyv[0]);
    for (int i=0; i<6; i++) {
      dyv[i] *= -3.0/2.0 / j2 * sp;
    }

    dev_vec_deriv_(dS, &dyv[6]);
    for (int i=6; i<12; i++) {
      dyv[i] *= -3.0/2.0 / j2 * sp;
    }

    double dS2[6];
    sub_vec(s, &alpha[6], 6, dS2);

    // q (0)
    dyv[12] = 0.0;

    // sa
    dyv[13] = -nT * pow(fabs(q), nT - 1.0) * copysign(1.0, q) / DT;
  }
  else {
    std::fill(dyv, dyv+nhist(), 0.0);
  }
  

  return 0;

}

// Flow rule
int YaguchiGr91FlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  std::fill(gv, gv+6, 0.0);

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  double Jn = J2_(dS);

  dev_vec(dS);

  if (Jn > 0.0) {
    for (int i=0; i<6; i++) {
      gv[i] = 3.0/2.0 * dS[i] / Jn;
    }
  }

  return 0;
}

int YaguchiGr91FlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  std::fill(dgv, dgv+36, 0.0);

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      if (i==j) {
        dgv[CINDEX(i,j,6)] = 2.0/3.0;
      }
      else {
        dgv[CINDEX(i,j,6)] = -1.0/3.0;
      }
    }
  }
  for (int i=3; i<6; i++) {
    dgv[CINDEX(i,i,6)] = 1.0;
  }

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS); 
  double j2 = J2_(dS);
  double mdS[6];
  dev_vec_deriv_(dS, mdS);
  dev_vec(dS);
  
  for (int i=0; i<6; i++) {
    dS[i] *= 3.0 / (2.0 * pow(j2,2.0));
  }

  outer_update_minus(dS, 6, mdS, 6, dgv);

  for (int i=0; i<36; i++) {
    dgv[i] *= 3.0/(2.0 * j2);
  }

  return 0;
}

int YaguchiGr91FlowRule::dg_da(const double * const s, const double * const alpha, double T,
             double * const dgv) const
{
  // Only the X terms have derivatives
  std::fill(dgv, dgv+(6*nhist()), 0.0);

  int nc = nhist();

  // Bizarrely this is the easiest way to do this
  double deriv[36];
  dg_ds(s, alpha, T, deriv);

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      dgv[CINDEX(i,(j+0),nc)] = -deriv[CINDEX(i,j,6)];
      dgv[CINDEX(i,(j+6),nc)] = -deriv[CINDEX(i,j,6)];
    }
  }

  return 0;
}

// Hardening rule
int YaguchiGr91FlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  double Jn = J2_(dS);
  
  double n[6];
  dev_vec(dS);
  for (int i=0; i<6; i++) {
    n[i] = 3.0/2.0 * dS[i] / Jn;
  }

  double C1i = C1(T);
  double a1i = a10(T) - alpha[12];

  double C2i = C2(T);
  double a2i = a2(T);

  // X1
  for (int i=0; i<6; i++) {
    hv[i+0] = C1i * (2.0/3.0 * a1i * n[i] - alpha[i+0]);
  }
  
  // X2
  for (int i=0; i<6; i++) {
    hv[i+6] = C2i * (2.0/3.0 * a2i * n[i] - alpha[i+6]);
  }

  // Q
  double qi = q(T);
  double di = d(T);
  hv[12] = di*(qi - alpha[12]);

  // sa
  double bri = br(T);
  double bhi = bh(T);
  double Ai = A(T);
  double Bi = B(T);
  double yi;
  y(s, alpha, T, yi);

  if (fabs(yi) > log_tol_) {
    double sas = Ai + Bi * log10(yi);
    
    if (sas < 0.0) {
      sas = 0.0;
    }
    double bi;
    if ((sas - alpha[13]) >= 0.0) {
      bi = bhi;
    }
    else {
      bi = bri;
    }
    hv[13] = bi * (sas - alpha[13]);
  }
  else {
    hv[13] = 0.0;
  }

  return 0;
}

int YaguchiGr91FlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  // Only the X terms have derivatives
  std::fill(dhv, dhv + nhist()*6, 0.0);

  double C1i = C1(T);
  double a1i = a10(T) - alpha[12];

  double C2i = C2(T);
  double a2i = a2(T);

  // Again, this is the easiest way to do this
  double deriv[36];
  dg_ds(s, alpha, T, deriv);

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((i+0),j,6)] = deriv[CINDEX(i,j,6)] * 2.0/3.0 * C1i * a1i;
      dhv[CINDEX((i+6),j,6)] = deriv[CINDEX(i,j,6)] * 2.0/3.0 * C2i * a2i;
    }
  }

  // The derivative of the rate wrt to the stress goes into the last row
  double bri = br(T);
  double bhi = bh(T);
  double Ai = A(T);
  double Bi = B(T);
  double yi;
  y(s, alpha, T, yi);
  if (fabs(yi) > log_tol_) {
    double sas = Ai + Bi * log10(yi);
    if (sas > 0.0) {
      double bi;
      if ((sas - alpha[13]) >= 0.0) {
        bi = bhi;
      }
      else {
        bi = bri;
      }
      dy_ds(s, alpha, T, &dhv[CINDEX(13,0,6)]);
      for (int i=0; i<6; i++) {
        dhv[CINDEX(13,i,6)] = bi * Bi / (yi * log(10.0)) * dhv[CINDEX(13,i,6)];
      }
    }
  }

  return 0;
}

int YaguchiGr91FlowRule::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  // Fair number of cross-terms are zero
  int nh = nhist();
  std::fill(dhv, dhv+(nh*nh), 0.0);

  // Generic X terms
  double deriv[6*nh];
  dg_da(s, alpha, T, deriv);
  double C1i = C1(T);
  double a1i = a10(T) - alpha[12];

  double C2i = C2(T);
  double a2i = a2(T);

  for (int i=0; i<6; i++) {
    for (int j=0; j<nh; j++) {
      dhv[CINDEX((i+0),j,nh)] = 2.0/3.0 * a1i * deriv[CINDEX(i,j,nh)];
      dhv[CINDEX((i+6),j,nh)] = 2.0/3.0 * a2i * deriv[CINDEX(i,j,nh)];
    }
  }

  for (int i=0; i<6; i++) {
    dhv[CINDEX((i+0),(i+0),nh)] -= 1.0;
    dhv[CINDEX((i+6),(i+6),nh)] -= 1.0;
  }


  for (int i=0; i<6; i++) {
    for (int j=0; j<nh; j++) {
      dhv[CINDEX((i+0),j,nh)] *= C1i;
      dhv[CINDEX((i+6),j,nh)] *= C2i;
    }
  }

  // X1 has an extra term for the time-varying a1
  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  double Jn = J2_(dS);
  
  double n[6];
  dev_vec(dS);
  for (int i=0; i<6; i++) {
    n[i] = 3.0/2.0 * dS[i] / Jn;
  }

  for (int i=0; i<6; i++) {
    dhv[CINDEX(i,12,nh)] -= C1i*2.0/3.0*n[i];
  }

  // Q is nice and easy
  double di = d(T);
  dhv[CINDEX(12,12,nh)] = -di;

  // There are two components to sa: the derivative of the rate wrt the history
  // and the derivative of the actual history term itself
  double bri = br(T);
  double bhi = bh(T);
  double Ai = A(T);
  double Bi = B(T);
  double yi;
  y(s, alpha, T, yi);
  
  if (fabs(yi) > log_tol_) {
    double sas = Ai + Bi * log10(yi);
    double bi;
    if ((sas - alpha[13]) >= 0.0) {
      bi = bhi;
    }
    else {
      bi = bri;
    }
    if (sas > 0.0) {
      dy_da(s, alpha, T, &dhv[CINDEX(13,0,nh)]);
      for (int i=0; i<nh; i++) {
        dhv[CINDEX(13,i,nh)] = bi * Bi / (yi * log(10.0)) * dhv[CINDEX(13,i,nh)];
      }
    }
    dhv[CINDEX(13,13,nh)] += -bi;
  }

  return 0;
}

// Hardening rule wrt to time
int YaguchiGr91FlowRule::h_time(const double * const s, 
                                const double * const alpha, double T,
                                double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
  
  double mi = m(T);

  // X1
  double g1i = g1(T);
  double J1 = J2_(&alpha[0]);
  for (int i=0; i<6; i++) {
    hv[i+0] = -g1i * pow(J1, mi-1.0) * alpha[i+0];
  }
  
  // X2
  double g2i = g2(T);
  double J2 = J2_(&alpha[6]);
  for (int i=0; i<6; i++) {
    hv[i+6] = -g2i * pow(J2, mi-1.0) * alpha[i+6];
  }

  return 0;
}

int YaguchiGr91FlowRule::dh_ds_time(const double * const s, 
                                    const double * const alpha, double T,
                                    double * const dhv) const
{
  // This is actually still zero
  std::fill(dhv, dhv+(nhist()*6), 0.0);
  return 0;
}

int YaguchiGr91FlowRule::dh_da_time(const double * const s, 
                                    const double * const alpha, double T,
                                    double * const dhv) const
{
  int nh = nhist();
  std::fill(dhv, dhv+(nh*nh), 0.0);

  // This is non-zero
  double mi = m(T);

  // X1
  double g1i = g1(T);
  double J1 = J2_(&alpha[0]);
  double X1[6];
  std::copy(&alpha[0], &alpha[6], X1);
  double X1d[6];
  dev_vec_deriv_(X1, X1d);
  dev_vec(X1);
  for (int i=0; i<6; i++) {
    X1[i] *= (mi-1.0) * pow(J1,mi-3.0) * 3.0 / 2.0;
  }

  double dX1[36];
  std::fill(dX1, dX1+36, 0.0);
  for (int i=0; i<6; i++) {
    dX1[CINDEX(i,i,6)] = pow(J1, mi-1.0);
  }
  outer_update(X1, 6, X1d, 6, dX1);
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((i+0),(j+0),nh)] = -g1i * dX1[CINDEX(i,j,6)];
    }
  }
  
  // X2
  double g2i = g2(T);
  double J2 = J2_(&alpha[6]);
  double X2[6];
  std::copy(&alpha[6], &alpha[12], X2);
  double X2d[6];
  dev_vec_deriv_(X2, X2d);
  dev_vec(X2);
  for (int i=0; i<6; i++) {
    X2[i] *= (mi-1.0) * pow(J2,mi-3.0) * 3.0 / 2.0;
  }

  double dX2[36];
  std::fill(dX2, dX2+36, 0.0);
  for (int i=0; i<6; i++) {
    dX2[CINDEX(i,i,6)] = pow(J2, mi-1.0);
  }
  outer_update(X2, 6, X2d, 6, dX2);
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((i+6),(j+6),nh)] = -g2i * dX2[CINDEX(i,j,6)];
    }
  }

  return 0;
}


// Properties...

double YaguchiGr91FlowRule::D(double T) const
{
  return pow(-1.119e-9 + 5.145e-11*T - 5.450e-14*T*T, -1.0/2.85);
}

double YaguchiGr91FlowRule::n(double T) const
{
  return 2.850;
}

double YaguchiGr91FlowRule::a10(double T) const
{
 return 2.082e3 - 8.110*T + 1.321e-2 * T*T - 7.278e-6 * T*T*T;
}

double YaguchiGr91FlowRule::C2(double T) const
{
  if (T < 798.0) {
    return 2.0e2;
  }
  else {
    return -2.992e3 + 4.000 * T;
  }
}

double YaguchiGr91FlowRule::a2(double T) const
{
  if (T < 773.0) {
    return 1.100e2;
  }
  else {
    return 3.373e3 - 7.622 * T + 4.400e-3 * T*T;
  }
}

double YaguchiGr91FlowRule::g1(double T) const
{
  if (T < 773.0) {
    return 5.064e-137 * exp(3.545e-1 * T);
  }
  else {
    return 5.336e-33 * exp(4.470e-2 * T);
  }
}

double YaguchiGr91FlowRule::g2(double T) const
{
  if (T < 773.0) {
    return 8.572e-108 * exp(2.771e-1 * T);
  }
  else {
    return 2.817e-11 - 7.538e-14 * T + 5.039e-17 * T*T;
  }
}

double YaguchiGr91FlowRule::m(double T) const
{
  if (T < 673.0) {
    return 1.200e1;
  }
  else if ((673.0 <= T) && (T < 773)) {
    return 5.766e2 * exp(-5.754e-3 * T);
  }
  else {
    return 6.750;
  }
}

double YaguchiGr91FlowRule::br(double T) const
{
  return 1.0e3;
}

double YaguchiGr91FlowRule::bh(double T) const
{
  return 1.065e3 - 1.404 * T + 4.417e-4 * T*T;
}

double YaguchiGr91FlowRule::A(double T) const
{
  if (T < 673.0) {
    return -1.841e2 + 2.075e-1 * T;
  }
  else if ((673.0 <= T) && (T < 823)) {
    return -4.799e2 + 1.262 * T - 9.133e-4 * T*T;
  }
  else {
    return -6.000e1;
  }
}

double YaguchiGr91FlowRule::B(double T) const
{
  if (T < 773.0) {
    return -6.340e1 + 6.850e-2 * T;
  }
  else {
    return -1.730e1;
  }
}

double YaguchiGr91FlowRule::d(double T) const
{
  return 3.208 - 1.010e-2 * T + 1.012e-5 * T*T;
}

double YaguchiGr91FlowRule::q(double T) const
{
  if (T < 673.0) {
    return 9.500e1;
  }
  else if ((673.0 <= T) && (T < 823.0)) {
    return -2.467e2 + 7.320e-1 * T - 3.333e-4 * T*T;
  }
  else {
    return 1.300e2;
  }
}

double YaguchiGr91FlowRule::C1(double T) const
{
  if (T < 673.0) {
    return 1.500e3;
  }
  else if ((673.0 <= T) && (T < 773.0)) {
    return -2.879e4 + 45.0 * T;
  }
  else {
    return 6.000e3;
  }
}

// Couple of helpers
double YaguchiGr91FlowRule::J2_(const double * const v) const
{
  double cv[6];
  std::copy(v, v+6, cv);
  dev_vec(cv);
  return sqrt(3.0/2.0 * dot_vec(cv, cv, 6));
}

void YaguchiGr91FlowRule::dev_vec_deriv_(const double * const a, 
                                        double * const b) const
{
  for (int i=0; i<3; i++) {
    b[i] = 2.0/3.0 * a[i];
    for (int j=0; j<3; j++) {
      if (i == j) continue;
      b[i] -= a[j]/3.0;
    }
  }
  for (int i=3; i<6; i++) {
    b[i] = a[i];
  }

}

} // namespace neml
