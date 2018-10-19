#include "models.h"

#include "nemlmath.h"
#include "nemlerror.h"

#include <cassert>
#include <limits>

namespace neml {

// NEMLModel_sd implementation
NEMLModel_sd::NEMLModel_sd(
    std::shared_ptr<LinearElasticModel> emodel,
    std::shared_ptr<Interpolate> alpha) :
      NEMLModel(), elastic_(emodel), alpha_(alpha)
{

}

NEMLModel_sd::~NEMLModel_sd()
{

}

size_t NEMLModel_sd::nstore() const
{
  return nhist(); // Get to this later...
}

int NEMLModel_sd::init_store(double * const store) const
{
  init_hist(store); // Get to this later

  return 0;
}

double NEMLModel_sd::alpha(double T) const
{
  return alpha_->value(T);
}

const std::shared_ptr<const LinearElasticModel> NEMLModel_sd::elastic() const
{
  return elastic_;
}

int NEMLModel_sd::elastic_strains(const double * const s_np1,
                                           double T_np1, 
                                           const double * const h_np1,
                                           double * const e_np1) const
{
  double S[36];
  elastic_->S(T_np1, S);
  mat_vec(S, 6, s_np1, 6, e_np1);

  return 0;
}

double NEMLModel_sd::bulk(double T) const
{
  return elastic_->K(T);
}

double NEMLModel_sd::shear(double T) const
{
  return elastic_->G(T);
}

int NEMLModel_sd::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  return 0;
}

// Implementation of small strain elasticity
SmallStrainElasticity::SmallStrainElasticity(
    std::shared_ptr<LinearElasticModel> elastic,
    std::shared_ptr<Interpolate> alpha) :
    NEMLModel_sd(elastic, alpha)
{

}

std::string SmallStrainElasticity::type()
{
  return "SmallStrainElasticity";
}

ParameterSet SmallStrainElasticity::parameters()
{
  ParameterSet pset(SmallStrainElasticity::type());

  pset.add_parameter<NEMLObject>("elastic");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          std::make_shared<ConstantInterpolate>(0.0));

  return pset;
}

std::unique_ptr<NEMLObject> SmallStrainElasticity::initialize(ParameterSet & params)
{
  return make_unique<SmallStrainElasticity>(
      params.get_object_parameter<LinearElasticModel>("elastic"),
      params.get_object_parameter<Interpolate>("alpha")
      ); 
}

size_t SmallStrainElasticity::nhist() const
{
  return 0;
}

int SmallStrainElasticity::init_hist(double * const hist) const
{
  return 0;
}

int SmallStrainElasticity::update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n)
{
  elastic_->C(T_np1, A_np1);
  mat_vec(A_np1, 6, e_np1, 6, s_np1);

  // Energy calculation (trapezoid rule)
  double de[6];
  double ds[6];
  sub_vec(e_np1, e_n, 6, de);
  add_vec(s_np1, s_n, 6, ds);
  for (int i=0; i<6; i++) ds[i] /= 2.0;
  u_np1 = u_n + dot_vec(ds, de, 6);
  p_np1 = p_n;

  return 0;
}

// Implementation of perfect plasticity
SmallStrainPerfectPlasticity::SmallStrainPerfectPlasticity(
    std::shared_ptr<LinearElasticModel> elastic,
    std::shared_ptr<YieldSurface> surface,
    std::shared_ptr<Interpolate> ys,
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter,
    bool verbose, int max_divide) :
      NEMLModel_sd(elastic, alpha),
      surface_(surface), ys_(ys),
      tol_(tol), miter_(miter), verbose_(verbose), max_divide_(max_divide)
{

}

std::string SmallStrainPerfectPlasticity::type()
{
  return "SmallStrainPerfectPlasticity";
}

ParameterSet SmallStrainPerfectPlasticity::parameters()
{
  ParameterSet pset(SmallStrainPerfectPlasticity::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("surface");
  pset.add_parameter<NEMLObject>("ys");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          std::make_shared<ConstantInterpolate>(0.0));
  pset.add_optional_parameter<double>("tol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<int>("max_divide", 8);

  return pset;
}

std::unique_ptr<NEMLObject> SmallStrainPerfectPlasticity::initialize(ParameterSet & params)
{
  return make_unique<SmallStrainPerfectPlasticity>(
      params.get_object_parameter<LinearElasticModel>("elastic"),
      params.get_object_parameter<YieldSurface>("surface"),
      params.get_object_parameter<Interpolate>("ys"),
      params.get_object_parameter<Interpolate>("alpha"),
      params.get_parameter<double>("tol"),
      params.get_parameter<int>("miter"),
      params.get_parameter<bool>("verbose"),
      params.get_parameter<int>("max_divide")
      ); 
}

int SmallStrainPerfectPlasticity::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Setup for substepping
  int nd = 0;                     // How many times we subdivided
  int tf = pow(2, max_divide_);   // Total integer step count
  int cm = tf;                    // Attempted step
  int cs = 0;                     // Integer step fraction

  double e_diff[6];
  for (int i=0; i<6; i++) e_diff[i] = e_np1[i] - e_n[i];
  double T_diff = T_np1 - T_n;
  double t_diff = t_np1 - t_n;

  double e_past[6];
  std::copy(e_n, e_n+6, e_past);
  // Ignore history, knowing it's blank
  double s_past[6];
  std::copy(s_n, s_n+6, s_past);
  double T_past = T_n;
  double t_past = t_n;
  double u_past = u_n;
  double p_past = p_n;

  double e_next[6];
  double s_next[6];
  double T_next;
  double t_next;
  double u_next;
  double p_next;

  while (cs < tf) {
    // Goal
    double sm = (double) (cs + cm) / (double) tf;
    for (int i=0; i<6; i++) e_next[i] = e_n[i] + sm * e_diff[i];
    T_next = T_n + sm * T_diff;
    t_next = t_n + sm * t_diff;

    int ier = update_substep_(e_next, e_past, T_next, T_past, t_next,
                              t_past, s_next, s_past, h_np1, h_n,
                              A_np1, u_next, u_past, p_next, p_past);

    // Subdivide
    if (ier != SUCCESS) {
      nd += 1;
      if (nd >= max_divide_) {
        return ier;
      }
      cm /= 2;
      continue;
    }

    // Next substep
    cs += cm;
    std::copy(e_next, e_next+6, e_past);
    std::copy(s_next, s_next+6, s_past);
    T_past = T_next;
    t_past = t_next;
    u_past = u_next;
    p_past = p_next;
  }

  // Final values
  std::copy(s_next, s_next+6, s_np1);
  u_np1 = u_next;
  p_np1 = p_next;

  return 0; 
}

int SmallStrainPerfectPlasticity::update_substep_(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Setup trial state
  SSPPTrialState ts;
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n, ts);

  // Check if this is an elastic state
  double fv;
  surface_->f(ts.s_tr, &ts.ys, T_np1, fv);
  if (fv < tol_) {
    std::copy(ts.s_tr, ts.s_tr+6, s_np1);
    std::copy(ts.C, ts.C+36, A_np1);

    p_np1 = p_n;
  }
  else {
    // Newton
    double * x = new double[nparams()];
    int ier = solve(this, x, &ts, tol_, miter_, verbose_);
    if (ier != SUCCESS) return ier;
    
    // Extract
    std::copy(x, x+6, s_np1);

    // Calculate tangent
    ier = calc_tangent_(ts, s_np1, x[6], A_np1);
    if (ier != SUCCESS) return ier;

    // Plastic work calculation
    double de[6];
    double ds[6];
    sub_vec(e_np1, e_n, 6, de);
    add_vec(s_np1, s_n, 6, ds);
    double ee_np1[6];
    mat_vec(ts.S, 6, s_np1, 6, ee_np1);
    sub_vec(de, ee_np1, 6, de);
    add_vec(de, ts.ee_n, 6, de);
    p_np1 = p_n + dot_vec(ds, de, 6) / 2.0;
    delete [] x;
  }

  // Energy calculation (trapezoid rule)
  double de[6];
  double ds[6];
  sub_vec(e_np1, e_n, 6, de);
  add_vec(s_np1, s_n, 6, ds);
  u_np1 = u_n + dot_vec(ds, de, 6) / 2.0;

  return 0;
}

size_t SmallStrainPerfectPlasticity::nhist() const
{
  return 0;
}

int SmallStrainPerfectPlasticity::init_hist(double * const hist) const
{
  return 0;
}

size_t SmallStrainPerfectPlasticity::nparams() const
{
  return 7;
}

int SmallStrainPerfectPlasticity::init_x(double * const x, TrialState * ts)
{
  SSPPTrialState * tss = static_cast<SSPPTrialState *>(ts);
  std::copy(tss->s_tr, tss->s_tr+6, x);
  x[6] = 0.0;

  return 0;
}

int SmallStrainPerfectPlasticity::RJ(
    const double * const x, TrialState * ts, double * const R,
    double * const J)
{
  SSPPTrialState * tss = static_cast<SSPPTrialState *>(ts);
  const double * const s_np1 = x;
  double dg = x[6];

  // Common things
  double fv;
  surface_->f(s_np1, &tss->ys, tss->T, fv);

  double df[6];
  surface_->df_ds(s_np1, &tss->ys, tss->T, df);

  double ddf[36];
  surface_->df_dsds(s_np1, &tss->ys, tss->T, ddf);

  // R1
  mat_vec(tss->S, 6, s_np1, 6, R);
  sub_vec(R, tss->e_np1, 6, R);
  add_vec(R, tss->e_n, 6, R);
  sub_vec(R, tss->ee_n, 6, R);
  for (int i=0; i<6; i++) R[i] += df[i] * dg;
  
  // R2
  R[6] = fv;

  // J11
  for (int i=0; i<36; i++) ddf[i] = ddf[i] * dg + tss->S[i];
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,7)] = ddf[CINDEX(i,j,6)];
    }
  }

  // J12
  for (int i=0; i<6; i++) J[CINDEX(i,6,7)] = df[i];

  // J21
  for (int i=0; i<6; i++) J[CINDEX(6,i,7)] = df[i];

  // J22
  J[CINDEX(6,6,7)] = 0.0;

  return 0;
}

// Getter
double SmallStrainPerfectPlasticity::ys(double T) const {
  return ys_->value(T);
}


// Make this public for ease of testing
int SmallStrainPerfectPlasticity::make_trial_state(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const double * const s_n, const double * const h_n,
    SSPPTrialState & ts)
{
  ts.ys = -ys_->value(T_np1);

  elastic_->S(T_np1, ts.S);

  std::copy(s_n, s_n+6, ts.s_n);

  double S_n[36];
  elastic_->S(T_n, S_n);
  mat_vec(S_n, 6, s_n, 6, ts.ee_n);

  double temp[6];
  elastic_->C(T_np1, ts.C);
  sub_vec(e_np1, e_n, 6, temp);

  add_vec(temp, ts.ee_n, 6, temp);
  mat_vec(ts.C, 6, temp, 6, ts.s_tr);

  ts.T = T_np1;

  std::copy(e_np1, e_np1+6, ts.e_np1);
  std::copy(e_n, e_n+6, ts.e_n);

  return 0;
}

int SmallStrainPerfectPlasticity::calc_tangent_(SSPPTrialState ts, 
                                                const double * const s_np1, 
                                                double dg, 
                                                double * const A_np1)
{
  // Useful
  double df[6];
  surface_->df_ds(s_np1, &ts.ys, ts.T, df);
  surface_->df_dsds(s_np1, &ts.ys, ts.T, A_np1);

  // Tangent calc
  for (int i=0; i<36; i++) A_np1[i] = ts.S[i] + dg * A_np1[i];
  invert_mat(A_np1, 6);

  double Bv[6];
  mat_vec(A_np1, 6, df, 6, Bv);

  double fact = dot_vec(Bv, df, 6);

  double Bvdiv[6];
  std::copy(Bv, Bv+6, Bvdiv);
  for (int i=0; i<6; i++) Bvdiv[i] /= fact;

  outer_update_minus(Bv, 6, Bvdiv, 6, A_np1);

  return 0;
}



// Implementation of small strain rate independent plasticity
//
SmallStrainRateIndependentPlasticity::SmallStrainRateIndependentPlasticity(
    std::shared_ptr<LinearElasticModel> elastic,
    std::shared_ptr<RateIndependentFlowRule> flow, 
    std::shared_ptr<Interpolate> alpha, double tol,
    int miter, bool verbose, double kttol, bool check_kt) :
      NEMLModel_sd(elastic, alpha),
      flow_(flow), tol_(tol), kttol_(kttol), miter_(miter),
      verbose_(verbose), check_kt_(check_kt)
{

}

std::string SmallStrainRateIndependentPlasticity::type()
{
  return "SmallStrainRateIndependentPlasticity";
}

ParameterSet SmallStrainRateIndependentPlasticity::parameters()
{
  ParameterSet pset(SmallStrainRateIndependentPlasticity::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("flow");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          std::make_shared<ConstantInterpolate>(0.0));
  pset.add_optional_parameter<double>("tol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<double>("kttol", 1.0e-2);
  pset.add_optional_parameter<bool>("check_kt", false);

  return pset;
}

std::unique_ptr<NEMLObject> SmallStrainRateIndependentPlasticity::initialize(ParameterSet & params)
{
  return make_unique<SmallStrainRateIndependentPlasticity>(
      params.get_object_parameter<LinearElasticModel>("elastic"),
      params.get_object_parameter<RateIndependentFlowRule>("flow"),
      params.get_object_parameter<Interpolate>("alpha"),
      params.get_parameter<double>("tol"),
      params.get_parameter<int>("miter"),
      params.get_parameter<bool>("verbose"),
      params.get_parameter<double>("kttol"),
      params.get_parameter<bool>("check_kt")
      ); 
}

size_t SmallStrainRateIndependentPlasticity::nhist() const
{
  return flow_->nhist();
}

int SmallStrainRateIndependentPlasticity::init_hist(double * const hist) const
{
  return flow_->init_hist(hist);
}

int SmallStrainRateIndependentPlasticity::update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n)
{
  // Setup and store the trial state for the solver
  SSRIPTrialState ts;
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n, ts);

  // Check to see if this is an elastic state
  double fv;
  flow_->f(ts.s_tr, &ts.h_tr[0], T_np1, fv);
  double dg;

  // If elastic, copy over and return
  if (fv < tol_) {
    std::copy(ts.s_tr, ts.s_tr+6, s_np1);
    std::copy(&ts.h_tr[0], &ts.h_tr[0]+flow_->nhist(), h_np1);
    std::copy(ts.C, ts.C+36, A_np1);
    dg = 0.0;

    p_np1 = p_n;
  }
  // Else solve and extract updated parameters from the solver vector
  else {
    double * x = new double[nparams()];
    int ier = solve(this, x, &ts, tol_, miter_, verbose_);
    if (ier != SUCCESS) return ier;

    // Extract solved parameters
    std::copy(x+6, x+6+flow_->nhist(), h_np1); // history
    dg = x[6+flow_->nhist()];
    double ee[6];
    sub_vec(e_np1, x, 6, ee);
    mat_vec(ts.C, 6, ee, 6, s_np1);

    // Complicated tangent calc...
    calc_tangent_(x, &ts, s_np1, h_np1, dg, A_np1);

    // Plastic work calculation
    double dep[6];
    double ds[6];
    add_vec(s_np1, s_n, 6, ds);
    sub_vec(x, ts.ep_tr, 6, dep);
    p_np1 = p_n + dot_vec(ds, dep, 6) / 2.0;
    
    delete [] x;
  }

  // Energy calculation (trapezoid rule)
  double de[6];
  double ds[6];
  sub_vec(e_np1, e_n, 6, de);
  add_vec(s_np1, s_n, 6, ds);
  u_np1 = u_n + dot_vec(ds, de, 6) / 2.0;
  
  // Check K-T and return
  return check_K_T_(s_np1, h_np1, T_np1, dg);

}

size_t SmallStrainRateIndependentPlasticity::nparams() const
{
  // My convention: plastic strain, history, consistency parameter
  return 6 + flow_->nhist() + 1;
}

int SmallStrainRateIndependentPlasticity::init_x(double * const x, TrialState * ts)
{
  SSRIPTrialState * tss = static_cast<SSRIPTrialState *>(ts);
  std::copy(tss->ep_tr, tss->ep_tr+6, x);
  std::copy(&tss->h_tr[0], &tss->h_tr[0]+flow_->nhist(), &x[6]);
  x[6+flow_->nhist()] = 0.0; // consistency parameter
  return 0;
}

int SmallStrainRateIndependentPlasticity::RJ(const double * const x, 
                                             TrialState * ts, 
                                             double * const R, double * const J)
{
  SSRIPTrialState * tss = static_cast<SSRIPTrialState *>(ts);

  // Setup from current state
  const double * const ep = &x[0];
  const double * const alpha  = &x[6];
  const double & dg = x[6+flow_->nhist()];
  double ee[6];
  sub_vec(tss->e_np1, ep, 6, ee);
  double s[6];
  mat_vec(tss->C, 6, ee, 6, s);

  // Residual calculation
  double g[6];
  flow_->g(s, alpha, tss->T, g); 
  double * h = new double[flow_->nhist()];
  flow_->h(s, alpha, tss->T, h);
  double f;
  flow_->f(s, alpha, tss->T, f);
  for (int i=0; i<6; i++) 
  {
    R[i] = -ep[i] + tss->ep_tr[i] + g[i] * dg;
  }
  for (size_t i=0; i<flow_->nhist(); i++) {
    R[i+6] = -alpha[i] + tss->h_tr[i] + h[i] * dg;
  }
  R[6+flow_->nhist()] = f;

  // Now the jacobian calculation...
  int n = nparams();
  int nh = flow_->nhist();
  
  // J11
  double gs[36];
  double J11[36];
  flow_->dg_ds(s, alpha, tss->T, gs);
  mat_mat(6, 6, 6, gs, tss->C, J11);
  for (int i=0; i<36; i++) J11[i] = -J11[i] * dg;
  for (int i=0; i<6; i++) J11[CINDEX(i,i,6)] -= 1.0;
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,n)] = J11[CINDEX(i,j,6)];
    }
  }
  
  // J12
  double * J12 = new double[6*nh];
  flow_->dg_da(s, alpha, tss->T, J12);
  for (int i=0; i<6*nh; i++) J12[i] *= dg;
  for (int i=0; i<6; i++) {
    for (int j=0; j < nh; j++) {
      J[CINDEX(i,(j+6),n)] = J12[CINDEX(i,j,nh)];
    }
  }
  delete [] J12;

  // J13
  for (int i=0; i<6; i++) {
    J[CINDEX(i,(6+nh),n)] = g[i];
  }

  // J21
  double * ha = new double[nh*6];
  double * J21 = new double[nh*6];
  flow_->dh_ds(s, alpha, tss->T, ha);
  mat_mat(nh, 6, 6, ha, tss->C, J21);
  for (int i=0; i<nh*6; i++) J21[i] = -J21[i] * dg;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,n)] = J21[CINDEX(i,j,6)];
    }
  }
  delete [] ha;
  delete [] J21;

  // J22
  double * J22 = new double[nh*nh];
  flow_->dh_da(s, alpha, tss->T, J22);
  for (int i=0; i<nh*nh; i++) J22[i] *= dg;
  for (int i=0; i<nh; i++) J22[CINDEX(i,i,nh)] -= 1.0;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<nh; j++) {
      J[CINDEX((i+6), (j+6), n)] = J22[CINDEX(i,j,nh)];
    }
  }
  delete [] J22;

  // J23
  for (int i=0; i<nh; i++) {
    J[CINDEX((i+6), (6+nh), n)] = h[i];
  }

  // J31
  double fs[6];
  double J31[6];
  flow_->df_ds(s, alpha, tss->T, fs);
  mat_vec_trans(tss->C, 6, fs, 6, J31);
  for (int i=0; i<6; i++) J31[i] = -J31[i];
  for (int i=0; i<6; i++) {
    J[CINDEX((6+nh), i, n)] = J31[i];
  }

  // J32
  double * J32 = new double[nh];
  flow_->df_da(s, alpha, tss->T, J32);
  for (int i=0; i<nh; i++) {
    J[CINDEX((6+nh), (i+6), n)] = J32[i];
  }
  delete [] J32;

  // J33
  J[CINDEX((6+nh), (6+nh), n)] = 0.0;
  
  delete [] h;

  return 0;
}

const std::shared_ptr<const LinearElasticModel> SmallStrainRateIndependentPlasticity::elastic() const
{
  return elastic_;
}

int SmallStrainRateIndependentPlasticity::make_trial_state(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const double * const s_n, const double * const h_n,
    SSRIPTrialState & ts)
{
  // Save e_np1
  std::copy(e_np1, e_np1+6, ts.e_np1);
  // ep_tr = ep_n
  double S_n[36];
  double ee_n[6];
  elastic_->S(T_n, S_n);
  mat_vec(S_n, 6, s_n, 6, ee_n);
  sub_vec(e_n, ee_n, 6, ts.ep_tr);
  // h_tr = h_n
  ts.h_tr.resize(flow_->nhist());
  std::copy(h_n, h_n+nhist(), ts.h_tr.begin());
  // Calculate the trial stress
  double ee[6];
  sub_vec(e_np1, ts.ep_tr, 6, ee);
  elastic_->C(T_np1, ts.C);
  mat_vec(ts.C, 6, ee, 6, ts.s_tr);
  // Store temp
  ts.T = T_np1;
  return 0;
}


int SmallStrainRateIndependentPlasticity::calc_tangent_(
    const double * const x, TrialState * ts, const double * const s_np1,
    const double * const h_np1, double dg, double * const A_np1)
{
  SSRIPTrialState * tss = static_cast<SSRIPTrialState *>(ts);

  double * R = new double[nparams()];
  double * J = new double[nparams()*nparams()];
  
  RJ(x, ts, R, J);

  delete [] R;

  int n = nparams();
  int nk = 6;
  int ne = nparams() - nk;
  
  double * Jkk = new double[nk*nk];
  for (int i=0; i<nk; i++) {
    for (int j=0; j<nk; j++) {
      Jkk[CINDEX(i,j,nk)] = J[CINDEX(i,j,n)];
    }
  }

  double * Jke = new double[nk*ne];
  for (int i=0; i<nk; i++) {
    for (int j=0; j<ne; j++) {
      Jke[CINDEX(i,j,ne)] = J[CINDEX(i,(j+nk),n)];
    }
  }

  double * Jek = new double[ne*nk];
  for (int i=0; i<ne; i++) {
    for (int j=0; j<nk; j++) {
      Jek[CINDEX(i,j,nk)] = J[CINDEX((i+nk),(j),n)];
    }
  }

  double * Jee = new double[ne*ne];
  for (int i=0; i<ne; i++) {
    for (int j=0; j<ne; j++) {
      Jee[CINDEX(i,j,ne)] = J[CINDEX((i+nk),(j+nk),n)];
    }
  }

  delete [] J;

  invert_mat(Jee, ne);
  
  double * A = new double[nk*6];
  double * B = new double[ne*6];
  
  int nh = flow_->nhist();

  double * dg_ds = new double[6*6];
  flow_->dg_ds(s_np1, h_np1, tss->T, dg_ds);
  double * dh_ds = new double[nh*6];
  flow_->dh_ds(s_np1, h_np1, tss->T, dh_ds);
  double df_ds[6];
  flow_->df_ds(s_np1, h_np1, tss->T, df_ds);

  mat_mat(6, 6, 6, dg_ds, tss->C, A);
  for (int i=0; i<nk*6; i++) A[i] *= dg;
  
  mat_mat(nh, 6, 6, dh_ds, tss->C, B);
  for (int i=0; i<nh*6; i++) B[i] *= dg;
  mat_vec_trans(tss->C, 6, df_ds, 6, &B[nh*6]);

  delete [] dg_ds;
  delete [] dh_ds;
  
  double * T1 = new double[ne*nk];
  mat_mat(ne, nk, ne, Jee, Jek, T1);
  double * T2 = new double[nk*nk];
  mat_mat(nk, nk, ne, Jke, T1, T2);
  for (int i=0; i<nk*nk; i++) T2[i] = Jkk[i] - T2[i];
  invert_mat(T2, nk);

  delete [] T1;

  double * T3 = new double[ne*nk];
  mat_mat(ne, nk, ne, Jee, B, T3);
  double * T4 = new double[nk*nk];
  mat_mat(nk, nk, ne, Jke, T3, T4);
  for (int i=0; i<nk*nk; i++) T4[i] -= A[i];

  delete [] T3;

  double dep[36];
  mat_mat(6, 6, 6, T2, T4, dep);

  delete [] T2;
  delete [] T4;

  for (int i=0; i<36; i++) dep[i] = -dep[i];
  for (int i=0; i<6; i++) dep[CINDEX(i,i,6)] += 1.0;

  mat_mat(6, 6, 6, tss->C, dep, A_np1);

  delete [] A;
  delete [] B;

  delete [] Jkk;
  delete [] Jke;
  delete [] Jek;
  delete [] Jee;

  return 0;
}

int SmallStrainRateIndependentPlasticity::check_K_T_(
    const double * const s_np1, const double * const h_np1, double T_np1, double dg)
{
  if (not check_kt_) {
    return 0;
  }

  double fv;
  flow_->f(s_np1, h_np1, T_np1, fv);

  if (fv > kttol_) {
    return KT_VIOLATION;
  }
  if (dg < -kttol_) {
    return KT_VIOLATION;
  }
  if ((dg > kttol_) && (fv > kttol_)) {
    return KT_VIOLATION;
  }

  return 0;
}

// Implement creep + plasticity
// Implementation of small strain rate independent plasticity
//
SmallStrainCreepPlasticity::SmallStrainCreepPlasticity(
    std::shared_ptr<LinearElasticModel> elastic,
    std::shared_ptr<NEMLModel_sd> plastic,
    std::shared_ptr<CreepModel> creep,
    std::shared_ptr<Interpolate> alpha, double tol,
    int miter, bool verbose, double sf) :
      NEMLModel_sd(elastic, alpha),
      plastic_(plastic), creep_(creep), tol_(tol), sf_(sf),
      miter_(miter), verbose_(verbose)
{

}

std::string SmallStrainCreepPlasticity::type()
{
  return "SmallStrainCreepPlasticity";
}

ParameterSet SmallStrainCreepPlasticity::parameters()
{
  ParameterSet pset(SmallStrainCreepPlasticity::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("plastic");
  pset.add_parameter<NEMLObject>("creep");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          std::make_shared<ConstantInterpolate>(0.0));
  pset.add_optional_parameter<double>("tol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<double>("sf", 1.0e6);

  return pset;
}

std::unique_ptr<NEMLObject> SmallStrainCreepPlasticity::initialize(ParameterSet & params)
{
  return make_unique<SmallStrainCreepPlasticity>(
      params.get_object_parameter<LinearElasticModel>("elastic"),
      params.get_object_parameter<NEMLModel_sd>("plastic"),
      params.get_object_parameter<CreepModel>("creep"),
      params.get_object_parameter<Interpolate>("alpha"),
      params.get_parameter<double>("tol"),
      params.get_parameter<int>("miter"),
      params.get_parameter<bool>("verbose"),
      params.get_parameter<double>("sf")
      ); 
}

size_t SmallStrainCreepPlasticity::nhist() const
{
  // The elastic-plastic strain + the plastic model history
  return plastic_->nhist() + 6;
}

int SmallStrainCreepPlasticity::init_hist(double * const hist) const
{
  std::fill(hist, hist+6, 0.0);
  return plastic_->init_hist(&hist[6]);
}

int SmallStrainCreepPlasticity::update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n)
{

  // Solve the system to get the update
  SSCPTrialState ts;
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n, ts);
  double * x = new double[nparams()];
  int ier = solve(this, x, &ts, tol_, miter_, verbose_);
  if (ier != 0) return ier;

  // Store the ep strain
  std::copy(x, x+6, h_np1);

  // Do the plastic update to get the new history and stress
  double A[36];
  ier =  plastic_->update_sd(x, ts.ep_strain, T_np1, T_n,
                             t_np1, t_n, s_np1, s_n,
                             &h_np1[6], &h_n[6],
                             A, u_np1, u_n, p_np1, p_n);
  if (ier != 0) return ier;
  delete [] x;

  // Do the creep update to get a tangent component
  double creep_old[6];
  double creep_new[6];
  double B[36];
  for (int i=0; i<6; i++) {
    creep_old[i] = e_n[i] - ts.ep_strain[i];
  }
  ier = creep_->update(s_np1, creep_new, creep_old, T_np1, T_n,
                 t_np1, t_n, B);
  if (ier != 0) return ier;

  // Form the relatively simple tangent
  ier = form_tangent_(A, B, A_np1);
  if (ier != 0) return ier;

  // Energy calculation (trapezoid rule)
  double de[6];
  double ds[6];
  sub_vec(e_np1, e_n, 6, de);
  add_vec(s_np1, s_n, 6, ds);
  u_np1 = u_n + dot_vec(ds, de, 6) / 2.0;
  
  // Extra dissipation from the creep material
  double dec[6];
  sub_vec(creep_new, creep_old, 6, dec);
  p_np1 += dot_vec(ds, dec, 6) / 2.0;
  
  return 0;
}

size_t SmallStrainCreepPlasticity::nparams() const
{
  // Just the elastic-plastic strain
  return 6;
}

int SmallStrainCreepPlasticity::init_x(double * const x, TrialState * ts)
{
  SSCPTrialState * tss = static_cast<SSCPTrialState*>(ts);

  // Start out at last step's value
  std::copy(tss->ep_strain, tss->ep_strain + 6, x);
  
  return 0;
}

int SmallStrainCreepPlasticity::RJ(const double * const x, TrialState * ts, 
                                   double * const R, double * const J)
{
  SSCPTrialState * tss = static_cast<SSCPTrialState*>(ts);

  int ier;

  // First update the elastic-plastic model
  double s_np1[6];
  double A_np1[36];
  std::vector<double> h_np1;
  h_np1.resize(plastic_->nhist());
  double u_np1, u_n;
  double p_np1, p_n;
  u_n = 0.0;
  p_n = 0.0;

  ier = plastic_->update_sd(x, tss->ep_strain, tss->T_np1, tss->T_n,
                      tss->t_np1, tss->t_n, s_np1, tss->s_n,
                      &h_np1[0], &(tss->h_n[0]), A_np1,
                      u_np1, u_n, p_np1, p_n);
  if (ier != 0) return ier;

  // Then update the creep strain
  double creep_old[6];
  double creep_new[6];
  double B[36];
  for (int i=0; i<6; i++) {
    creep_old[i] = tss->e_n[i] - tss->ep_strain[i];
  }
  ier = creep_->update(s_np1, creep_new, creep_old, tss->T_np1, tss->T_n,
                 tss->t_np1, tss->t_n, B);
  if (ier != 0) return ier;

  // Form the residual
  for (int i=0; i<6; i++) {
    R[i] = (x[i] + creep_new[i] - tss->e_np1[i]) * sf_;
  }
  
  // The Jacobian is a straightforward combination of the two derivatives
  ier = mat_mat(6, 6, 6, B, A_np1, J);
  for (int i=0; i<6; i++) J[CINDEX(i,i,6)] += 1.0;
  for (int i=0; i<36; i++) J[i] *= sf_;

  return ier;
}

int SmallStrainCreepPlasticity::make_trial_state(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const double * const s_n, const double * const h_n,
    SSCPTrialState & ts)
{
  int nh = plastic_->nhist();
  ts.h_n.resize(nh);

  std::copy(e_np1, e_np1+6, ts.e_np1);
  std::copy(e_n, e_n+6, ts.e_n);
  std::copy(s_n, s_n+6, ts.s_n);
  ts.T_n = T_n;
  ts.T_np1 = T_np1;
  ts.t_n = t_n;
  ts.t_np1 = t_np1;
  std::copy(h_n + 6, h_n + 6 + nh, ts.h_n.begin());

  std::copy(h_n, h_n+6, ts.ep_strain);

  return 0;
}

int SmallStrainCreepPlasticity::form_tangent_(
    double * const A, double * const B, double * const A_np1)
{
  // Okay, what we really want to do is
  // (A^-1 + B)^-1
  // BUT A can be singular (i.e. if it's perfectly plastic)
  // The suspicion is that the whole thing probably isn't singular
  // even if A is.
  //
  // It turns out there is a formula for the inverse of the sum of two
  // matrices in Henderson and Searle (1981):
  // (A+B)^-1 = A^-1 - A^-1 B (I + A^-1 B)^-1 A^-1
  // which does not require *B* to be nonsingular
  // This reduces to
  // A - A B (I + A B)^-1 A
  // for our case, but of course it's not necessarily correct as our
  // A may be singular.
  //
  // That said, it seems to work quite nicely.
  //
  double C[36];
  mat_mat(6,6,6,A,B,C);
  for (int i=0; i<6; i++) C[CINDEX(i,i,6)] += 1.0;
  invert_mat(C, 6);

  double D[36];
  mat_mat(6,6,6,C,A,D);
  mat_mat(6,6,6,B,D,C);
  mat_mat(6,6,6,A,C,D);

  std::copy(A,A+36,A_np1);
  for (int i=0; i<36; i++) A_np1[i] -= D[i];

  return 0;
}

int SmallStrainCreepPlasticity::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  plastic_->set_elastic_model(emodel);
  return 0;
}

// Start general integrator implementation
GeneralIntegrator::GeneralIntegrator(std::shared_ptr<LinearElasticModel> elastic,
                                     std::shared_ptr<GeneralFlowRule> rule,
                                     std::shared_ptr<Interpolate> alpha,
                                     double tol, int miter,
                                     bool verbose, int max_divide) :
    NEMLModel_sd(elastic, alpha),
    rule_(rule), tol_(tol), miter_(miter), max_divide_(max_divide),
    verbose_(verbose) 
{

}

std::string GeneralIntegrator::type()
{
  return "GeneralIntegrator";
}

ParameterSet GeneralIntegrator::parameters()
{
  ParameterSet pset(GeneralIntegrator::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("rule");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          std::make_shared<ConstantInterpolate>(0.0));
  pset.add_optional_parameter<double>("tol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<int>("max_divide", 8);

  return pset;
}

std::unique_ptr<NEMLObject> GeneralIntegrator::initialize(ParameterSet & params)
{
  return make_unique<GeneralIntegrator>(
      params.get_object_parameter<LinearElasticModel>("elastic"),
      params.get_object_parameter<GeneralFlowRule>("rule"),
      params.get_object_parameter<Interpolate>("alpha"),
      params.get_parameter<double>("tol"),
      params.get_parameter<int>("miter"),
      params.get_parameter<bool>("verbose"),
      params.get_parameter<int>("max_divide")
      ); 
}

int GeneralIntegrator::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Setup for substepping
  int nd = 0;                   // Number of times we divided
  int tf = pow(2,max_divide_);  // Total integer step, to avoid floating math
  int cm = tf;                  // Current attempted step
  int cs = 0;                   // Current integer proportion of step completed
  
  // Total differences over step
  double e_diff[6];
  for (int i=0; i<6; i++) e_diff[i] = e_np1[i] - e_n[i];
  double T_diff = T_np1 - T_n;
  double t_diff = t_np1 - t_n;

  // Previous values as we go along, (initialize to step n)
  double e_past[6];
  std::copy(e_n, e_n+6, e_past);
  double * const h_past = new double [nhist()];
  std::copy(h_n, h_n+nhist(), h_past);
  double s_past[6];
  std::copy(s_n, s_n+6, s_past);
  double T_past = T_n;
  double t_past = t_n;
  
  // Current goal as we go along
  double e_next[6];
  double * const h_next = new double [nhist()];
  double s_next[6];
  double T_next;
  double t_next;
  
  while (cs < tf) {
    // Figure out our float step multiplier
    double sm = (double) (cs + cm) / (double) tf;

    // Figure out our goals for this increment
    for (int i=0; i<6; i++) e_next[i] = e_n[i] + sm * e_diff[i];
    T_next = T_n + sm * T_diff;
    t_next = t_n + sm * t_diff;

    // Set trial state
    GITrialState ts;
    make_trial_state(e_next, e_past, T_next, T_past, t_next, t_past, 
                     s_past, h_past, ts);

    // Solve for x
    double * x = new double[nparams()];
    int ier = solve(this, x, &ts, tol_, miter_, verbose_);

    // Decide what to do if we fail
    if (ier != SUCCESS) {
      // Subdivide the step
      nd += 1;
      cm /= 2;
      if (verbose_) {
        std::cout << "Substepping:" << std::endl;
        std::cout << "New step fraction " << ((double) cm / (double) tf) << std::endl;
        std::cout << "New integer step " << cm << std::endl;
        std::cout << "Step integer count " << cs << "/" << tf << std::endl;
      }
      // Check if we exceeded our subdivision limit
      if (nd == max_divide_) {
        if (verbose_) {
          std::cout << "Substepping failed..." << std::endl;
        }
        delete [] h_past;
        delete [] h_next;
        return ier;
      }
      continue;
    }

    // Extract solved parameters
    std::copy(x, x+6, s_next);
    std::copy(x+6, x+6+nhist(), h_next);

    delete [] x;

    // Increment next step
    cs += cm;
    std::copy(e_next, e_next+6, e_past);
    std::copy(h_next, h_next+nhist(), h_past);
    std::copy(s_next, s_next+6, s_past);
    T_past = T_next;
    t_past = t_next;

  }

  // Extract final values
  std::copy(s_next, s_next+6, s_np1);
  std::copy(h_next, h_next+nhist(), h_np1);
  
  // Free memory
  delete [] h_past;
  delete [] h_next;

  // Get tangent over full step
  GITrialState ts;
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n, ts);
  double * y = new double[nparams()];
  std::copy(s_np1, s_np1+6, y);
  std::copy(h_np1, h_np1+nhist(), &y[6]);
  
  calc_tangent_(y, &ts, A_np1);

  delete [] y;
  
  // Energy calculation (trapezoid rule)
  double de[6];
  double ds[6];
  sub_vec(e_np1, e_n, 6, de);
  add_vec(s_np1, s_n, 6, ds);
  for (int i=0; i<6; i++) ds[i] /= 2.0;
  u_np1 = u_n + dot_vec(ds, de, 6);

  // Need a special call
  double p_dot_np1;
  rule_->work_rate(s_np1, h_np1, ts.e_dot, T_np1, ts.Tdot, p_dot_np1);
  double p_dot_n;
  rule_->work_rate(s_n, h_n, ts.e_dot, T_n, ts.Tdot, p_dot_n);
  p_np1 = p_n + (p_dot_np1 + p_dot_n)/2.0 * ts.dt;

  return 0;
}

size_t GeneralIntegrator::nhist() const
{
  return rule_->nhist();
}

int GeneralIntegrator::init_hist(double * const hist) const
{
  return rule_->init_hist(hist);
}

size_t GeneralIntegrator::nparams() const
{
  return 6 + nhist();
}

int GeneralIntegrator::init_x(double * const x, TrialState * ts)
{
  GITrialState * tss = static_cast<GITrialState*>(ts);
  std::copy(tss->s_n, tss->s_n+6, x);
  std::copy(tss->h_n.begin(), tss->h_n.end(), &x[6]);

  return 0;
}

int GeneralIntegrator::RJ(const double * const x, TrialState * ts,
                          double * const R, double * const J)
{
  GITrialState * tss = static_cast<GITrialState*>(ts);

  // Setup
  double s_mod[6];
  std::copy(x, x+6, s_mod);
  if (norm2_vec(x, 6) < std::numeric_limits<double>::epsilon()) {
    s_mod[0] = 2.0 * std::numeric_limits<double>::epsilon();
  }
  const double * const h_np1 = &x[6];
  
  // Helps with vectorization
  // Really as I declared both const this shouldn't be necessary but hey
  // I don't design optimizing compilers for a living
  int nhist = this->nhist();
  int nparams = this->nparams();

  // Residual calculation
  rule_->s(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, R);
  for (int i=0; i<6; i++) {
    R[i] = -s_mod[i] + tss->s_n[i] + R[i] * tss->dt;
  }
  rule_->a(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, &R[6]);
  for (int i=0; i<nhist; i++) {
    R[i+6] = -h_np1[i] + tss->h_n[i] + R[i+6] * tss->dt;
  }

  // Jacobian calculation
  double J11[36];
  rule_->ds_ds(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, J11);
  for (int i=0; i<36; i++) J11[i] *= tss->dt;
  for (int i=0; i<6; i++) J11[CINDEX(i,i,6)] -= 1.0;
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,nparams)] = J11[CINDEX(i,j,6)];
    }
  }

  double * J12 = new double[6*nhist];
  rule_->ds_da(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, J12);
  for (int i=0; i<6; i++) {
    for (int j=0; j<nhist; j++) {
      J[CINDEX(i,(j+6),nparams)] = J12[CINDEX(i,j,nhist)] * tss->dt;
    }
  }

  delete [] J12;

  double * J21 = new double[nhist*6];
  rule_->da_ds(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, J21);
  for (int i=0; i<nhist; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,nparams)] = J21[CINDEX(i,j,6)] * tss->dt;
    }
  }
  delete [] J21;

  double * J22 = new double[nhist*nhist];
  rule_->da_da(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, J22);

  // More vectorization
  double dt = tss->dt;
  for (int i=0; i<nhist*nhist; i++) J22[i] *= dt;
  for (int i=0; i<nhist; i++) J22[CINDEX(i,i,nhist)] -= 1.0;

  for (int i=0; i<nhist; i++) {
    for (int j=0; j<nhist; j++) {
      J[CINDEX((i+6),(j+6),nparams)] = J22[CINDEX(i,j,nhist)];
    }
  }
  
  delete [] J22;

  return 0;
}


int GeneralIntegrator::make_trial_state(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const double * const s_n, const double * const h_n,
    GITrialState & ts)
{
  // Basic
  ts.dt = t_np1 - t_n;
  ts.T = T_np1;

  // Rates, looking out for divide-by-zero
  if (ts.dt > 0.0) {
    ts.Tdot = (T_np1 - T_n) / ts.dt;
    for (int i=0; i<6; i++) {
      ts.e_dot[i] = (e_np1[i] - e_n[i]) / ts.dt;
    }
  }
  else {
    ts.Tdot = 0.0;
    std::fill(ts.e_dot, ts.e_dot+6, 0.0);
  }
  
  // Trial stress
  std::copy(s_n, s_n+6, ts.s_n);

  // Trial history
  ts.h_n.resize(nhist());
  std::copy(h_n, h_n+nhist(), ts.h_n.begin());

  return 0;
}

int GeneralIntegrator::calc_tangent_(const double * const x, TrialState * ts, 
                                     double * const A_np1)
{
  // Quick note: I'm leaving  out a few dts that cancel in the end -- 
  // no point in tempting fate for small time increments
 
  GITrialState * tss = static_cast<GITrialState*>(ts);

  // Setup
  double s_mod[36];
  std::copy(x, x+6, s_mod);
  if (norm2_vec(x, 6) < std::numeric_limits<double>::epsilon()) {
    s_mod[0] = 2.0 * std::numeric_limits<double>::epsilon();
  }
  const double * const h_np1 = &x[6];

  // Vectorization
  int nhist = this->nhist();

  // Call for extra derivatives
  double A[36];
  rule_->ds_de(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, A);
  double * B = new double[nhist*6];
  rule_->da_de(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, B);

  // Call for the jacobian
  double * R = new double[nparams()];
  double * J = new double[nparams()*nparams()];
  RJ(x, ts, R, J);
  delete [] R;

  // Separate blocks...
  int n = nparams();
  
  double * J11 = new double[6*6];
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J11[CINDEX(i,j,6)] = J[CINDEX(i,j,n)];
    }
  }

  double * J12 = new double[6*nhist];
  for (int i=0; i<6; i++) {
    for (int j=0; j<nhist; j++) {
      J12[CINDEX(i,j,nhist)] = J[CINDEX(i,(j+6),n)];
    }
  }

  double * J21 = new double[nhist*6];
  for (int i=0; i<nhist; i++) {
    for (int j=0; j<6; j++) {
      J21[CINDEX(i,j,6)] = J[CINDEX((i+6),(j),n)];
    }
  }

  double * J22 = new double[nhist*nhist];
  for (int i=0; i<nhist; i++) {
    for (int j=0; j<nhist; j++) {
      J22[CINDEX(i,j,nhist)] = J[CINDEX((i+6),(j+6),n)];
    }
  }

  delete [] J;

  // Only need the inverse
  invert_mat(J22, nhist);

  // Start multiplying through
  double * T1 = new double[nhist*6];
  mat_mat(nhist, 6, nhist, J22, J21, T1);
  double * T2 = new double[6*6];
  mat_mat(6, 6, nhist, J12, T1, T2);
  for (int i=0; i<6*6; i++) T2[i] = J11[i] - T2[i];
  invert_mat(T2, 6);

  delete [] T1;

  double * T3 = new double[nhist*6];
  mat_mat(nhist, 6, nhist, J22, B, T3);
  double * T4 = new double[6*6];
  mat_mat(6, 6, nhist, J12, T3, T4);
  for (int i=0; i<6*6; i++) T4[i] -= A[i];

  delete [] T3;

  mat_mat(6, 6, 6, T2, T4, A_np1);

  delete [] T2;
  delete [] T4;

  delete [] B;
  delete [] J11;
  delete [] J12;
  delete [] J21;
  delete [] J22;

  return 0;
}

int GeneralIntegrator::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  rule_->set_elastic_model(emodel);
  return 0;
}

// Start KMRegimeModel
KMRegimeModel::KMRegimeModel(std::shared_ptr<LinearElasticModel> emodel,
                             std::vector<std::shared_ptr<NEMLModel_sd>> models,
                             std::vector<double> gs,
                             double kboltz, double b, double eps0,
                             std::shared_ptr<Interpolate> alpha) :
    NEMLModel_sd(emodel, alpha), models_(models), gs_(gs),
    kboltz_(kboltz), b_(b), eps0_(eps0)
{

}

std::string KMRegimeModel::type()
{
  return "KMRegimeModel";
}

ParameterSet KMRegimeModel::parameters()
{
  ParameterSet pset(KMRegimeModel::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<std::vector<NEMLObject>>("models");
  pset.add_parameter<std::vector<double>>("gs");
  pset.add_parameter<double>("kboltz");
  pset.add_parameter<double>("b");
  pset.add_parameter<double>("eps0");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          std::make_shared<ConstantInterpolate>(0.0));

  return pset;
}

std::unique_ptr<NEMLObject> KMRegimeModel::initialize(ParameterSet & params)
{
  return make_unique<KMRegimeModel>(
      params.get_object_parameter<LinearElasticModel>("elastic"),
      params.get_object_parameter_vector<NEMLModel_sd>("models"),
      params.get_parameter<std::vector<double>>("gs"),
      params.get_parameter<double>("kboltz"),
      params.get_parameter<double>("b"),
      params.get_parameter<double>("eps0"),
      params.get_object_parameter<Interpolate>("alpha")
      ); 
}

int KMRegimeModel::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Calculate activation energy
  double g = activation_energy_(e_np1, e_n, T_np1, t_np1, t_n);

  // Note this relies on everything being sorted.  You probably want to
  // error check at some point
  for (size_t i=0; i<gs_.size(); i++) {
    if (g < gs_[i]) {
      return models_[i]->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n, 
                                   s_np1, s_n, h_np1, h_n, A_np1, u_np1, u_n,
                                   p_np1, p_n);
    }
  }
  return models_.back()->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n,
                                   s_np1, s_n, h_np1, h_n, A_np1, u_np1, u_n,
                                   p_np1, p_n);
}

size_t KMRegimeModel::nhist() const
{
  return models_[0]->nhist();
}

int KMRegimeModel::init_hist(double * const hist) const
{
  return models_[0]->init_hist(hist);
}

double KMRegimeModel::activation_energy_(const double * const e_np1, 
                                         const double * const e_n,
                                         double T_np1,
                                         double t_np1, double t_n)
{
  double dt = t_np1 - t_n;

  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  for (int i=0; i<6; i++) de[i] /= dt;
  double rate = sqrt(2.0/3.0) * norm2_vec(de, 6);
  double mu = elastic_->G(T_np1);
  
  return kboltz_ * T_np1 / (mu* pow(b_, 3.0)) * log(eps0_ / rate);
}

int KMRegimeModel::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    (*it)->set_elastic_model(emodel);
  }
  return 0;
}

} // namespace neml
