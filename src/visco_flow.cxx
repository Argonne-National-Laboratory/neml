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

ChabocheFlowRule::ChabocheFlowRule(std::shared_ptr<YieldSurface> surface,
                                   std::shared_ptr<NonAssociativeHardening> hardening,
                                   std::shared_ptr<FluidityModel> fluidity,
                                   double n) :
    surface_(surface), hardening_(hardening), fluidity_(fluidity), n_(n)
{
  
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
    double eta = sqrt(2.0/3.0) * fluidity_->eta(alpha[0]);
    yv = sqrt(3.0/2.0) * pow(fv/eta, n_);
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
    double eta = sqrt(2.0/3.0) * fluidity_->eta(alpha[0]);
    double mv = sqrt(3.0/2.0) * pow(fv/eta, n_ - 1.0) * n_ / eta;
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

    double eta = sqrt(2.0/3.0) * fluidity_->eta(alpha[0]);
    double mv = sqrt(3.0/2.0) * pow(fv/eta, n_ - 1.0) * n_ / eta;
    for (int i=0; i<nhist(); i++) dyv[i] *= mv;

    double mv2 = -sqrt(3.0/2.0) * fv * pow(fv/eta, n_ - 1.0) * n_ / (eta * eta);
    double deta = sqrt(2.0/3.0) * fluidity_->deta(alpha[0]);
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


YaguchiGr91FlowRule::YaguchiGr91FlowRule()
{
  
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

  return 0;
}

int YaguchiGr91FlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{

  return 0;

}

// Flow rule
int YaguchiGr91FlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  double Jn = J2_(dS);

  for (int i=0; i<6; i++) {
    gv[i] = 3.0/2.0 * dS[i] / Jn;
  }

  return 0;
}

int YaguchiGr91FlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{

  return 0;
}

int YaguchiGr91FlowRule::dg_da(const double * const s, const double * const alpha, double T,
             double * const dgv) const
{

  return 0;
}

// Hardening rule
int YaguchiGr91FlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  double Jn = J2_(dS);
  
  double n[6];
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
  double sas = Ai + Bi * log10(yi);
  if (sas > 0.0) {
    sas = fabs(sas);
  }
  else {
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

  return 0;
}

int YaguchiGr91FlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return 0;
}

int YaguchiGr91FlowRule::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return 0;
}

// Hardening rule wrt to time
int YaguchiGr91FlowRule::h_time(const double * const s, 
                                const double * const alpha, double T,
                                double * const hv) const
{
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
  return 0;
}

int YaguchiGr91FlowRule::dh_da_time(const double * const s, 
                                    const double * const alpha, double T,
                                    double * const dhv) const
{
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

} // namespace neml
