#include "models.h"

#include "math/nemlmath.h"
#include "nemlerror.h"

#include <cassert>
#include <limits>
#include <iostream>
#include <fstream>

namespace neml {

NEMLModel::NEMLModel(ParameterSet & params) :
    HistoryNEMLObject(params)
{

}

void NEMLModel::save(std::string file_name, std::string model_name)
{
  std::string representation = serialize(model_name, "materials");
  
  std::ofstream outfile(file_name);
  outfile << representation;
  outfile.close();
}

void NEMLModel::populate_hist(History & h) const
{
  populate_static(h);
  populate_state(h);
}

void NEMLModel::init_hist(History & h) const
{
  h.zero(); // This shuts up valgrind errors.
  init_static(h);
  init_state(h);
}

void NEMLModel::populate_static(History & h) const
{
  // Nothing by default
}

void NEMLModel::init_static(History & h) const
{
  // Nothing by default
}

double NEMLModel::get_damage(const double *const h_np1) 
{ 
  return 0.0; 
}

bool NEMLModel::should_del_element(const double *const h_np1)
{
  return false;
}

bool NEMLModel::is_damage_model() const 
{
  return false;
}

std::tuple<History,History> NEMLModel::split_state(const History & h) const
{
  History h1 = h.split(stored_static_.items());
  History h2 = h.split(stored_static_.items(), false);

  return std::tie(h1,h2);
}

size_t NEMLModel::nstate() const
{
  return stored_state_.size();
}

size_t NEMLModel::nstatic() const
{
  return stored_static_.size();
}

void NEMLModel::cache_history_()
{
  HistoryNEMLObject::cache_history_();
  populate_state(stored_state_);
  populate_static(stored_static_);
}

std::vector<std::string> NEMLModel::report_internal_variable_names() const
{
  return gather_blank_history_().formatted_names();
}

// NEMLModel_sd implementation
NEMLModel_sd::NEMLModel_sd(ParameterSet & params) :
      NEMLModel(params), 
      elastic_(params.get_object_parameter<LinearElasticModel>("elastic")), 
      alpha_(params.get_object_parameter<Interpolate>("alpha")),
      truesdell_(params.get_parameter<bool>("truesdell"))
{

}

void NEMLModel_sd::update_ld_inc(
    const double * const d_np1, const double * const d_n,
    const double * const w_np1, const double * const w_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1, double * const B_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  double base_A_np1[36];
  
  double D[6];
  double W[3];
  double dS[6];

  sub_vec(d_np1, d_n, 6, D);
  sub_vec(w_np1, w_n, 3, W);
  
  // If not the Truesdell rate then use the Jaumann rate
  if (not truesdell_) {
    std::fill(D, D+6, 0.0);
  }
  
  update_sd(d_np1, d_n, T_np1, T_n, t_np1, t_n, &h_np1[0], &h_n[0],
                   h_np1, h_n, base_A_np1, u_np1, u_n, p_np1, p_n);

  sub_vec(&h_np1[0], &h_n[0], 6, dS);  
 
  truesdell_update_sym(D, W, s_n, dS, s_np1);

  calc_tangent_(D, W, base_A_np1, s_np1, A_np1, B_np1);
}

void NEMLModel_sd::update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n)
{
  update_sd_actual(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n,
            &h_np1[nstatic()], &h_n[nstatic()], 
            A_np1, u_np1, u_n, p_np1, p_n);
}

void NEMLModel_sd::init_static(History & h) const
{
  h.get<Symmetric>(prefix("small_stress")) = Symmetric::zero();
}

void NEMLModel_sd::populate_static(History & h) const
{
  // This is a legacy from the old flat vector system.  Composite
  // expect to only have this added one time.
  if (!h.contains("small_stress"))
    h.add<Symmetric>(prefix("small_stress"));
}

double NEMLModel_sd::alpha(double T) const
{
  return alpha_->value(T);
}

const std::shared_ptr<const LinearElasticModel> NEMLModel_sd::elastic() const
{
  return elastic_;
}

void NEMLModel_sd::elastic_strains(const double * const s_np1,
                                           double T_np1, 
                                           const double * const h_np1,
                                           double * const e_np1) const
{
  double S[36];
  elastic_->S(T_np1, S);
  mat_vec(S, 6, s_np1, 6, e_np1);

}

void NEMLModel_sd::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
}

void NEMLModel_sd::calc_tangent_(const double * const D, const double * const W,
                                const double * const C, const double * const S,
                                double * const A, double * const B)
{
  double J[81];
  double O[81];

  truesdell_mat(D, W, J);
  invert_mat(J, 9);

  truesdell_tangent_outer(S, O);

  double F[81];
  mandel2full(C, F);
  double dL[81];
  mat_mat(9,9,9, F, idsym, dL);

  double T1[81];
  for (int i=0; i<81; i++) T1[i] = dL[i] + O[i];

  double L[81];
  mat_mat(9,9,9, J, T1, L);

  double Dpart[81];
  double Wpart[81];

  mat_mat(9,9,9, L, idsym, Dpart);
  mat_mat(9,9,9, L, idskew, Wpart);

  full2mandel(Dpart, A);
  full2skew(Wpart, B);

  // IMPORTANT TODO: go back and find where you dropped the factor of 2...
  for (int i=0; i<18; i++) B[i] *= 2.0;

}

// NEMLModel_ldi implementation
NEMLModel_ldi::NEMLModel_ldi(ParameterSet & params) :
      NEMLModel(params)
{

}

void NEMLModel_ldi::update_sd(
   const double * const e_np1, const double * const e_n,
   double T_np1, double T_n,
   double t_np1, double t_n,
   double * const s_np1, const double * const s_n,
   double * const h_np1, const double * const h_n,
   double * const A_np1,
   double & u_np1, double u_n,
   double & p_np1, double p_n)
{
  double W[3] = {0,0,0};
  double B[18];
  update_ld_inc(e_np1, e_n, W, W, T_np1, T_n, t_np1, t_n,
                       s_np1, s_n, h_np1, h_n, A_np1, B, u_np1, u_n,
                       p_np1, p_n);
}

SubstepModel_sd::SubstepModel_sd(ParameterSet & params) :
    NEMLModel_sd(params), 
    rtol_(params.get_parameter<double>("rtol")),
    atol_(params.get_parameter<double>("atol")),
    miter_(params.get_parameter<int>("miter")),
    verbose_(params.get_parameter<bool>("verbose")), 
    linesearch_(params.get_parameter<bool>("linesearch")),
    max_divide_(params.get_parameter<int>("max_divide")),
    force_divide_(params.get_parameter<bool>("force_divide"))
{

}

void SubstepModel_sd::update_sd_actual(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
// Setup the substep parameters
  int nd = 0;                     // Number of times we subdivided
  int tf = pow(2, max_divide_);   // Total integer step count
  int cm = tf;                    // Attempted integer step (sf = cm/tf)
  int cs = 0;                     // Accumulated steps
  
  // Strain, time, and temperature increments
  double e_diff[6];
  sub_vec(e_np1, e_n, 6, e_diff);
  double T_diff = T_np1 - T_n;
  double t_diff = t_np1 - t_n;

  // Previous subincrement quantities
  double e_past[6];
  std::copy(e_n, e_n+6, e_past);
  double s_past[6];
  std::copy(s_n, s_n+6, s_past);
  double * h_past = new double [nstate()];
  std::copy(h_n, h_n+nstate(), h_past);
  double T_past = T_n;
  double t_past = t_n;
  double u_past = u_n;
  double p_past = p_n;

  // Targets
  double e_next[6];
  double T_next;
  double t_next;
  
  // Storage for the local A matrix
  double * A_inc = new double[nparams() * nparams()];
  double * A_old = new double[nparams() * 6];
  double * A_new = new double[nparams() * 6];
  double * E_inc = new double[nparams() * 6];

  std::fill(A_old, A_old+(nparams()*6), 0.0);

  while (cs < tf) {
    // targets
    double sm = (double) (cs + cm) / (double) tf;
    double sf = (double) cm / (double) tf;
    for (size_t i = 0; i<6; i++) e_next[i] = e_n[i] + sm * e_diff[i];
    T_next = T_n + sm * T_diff;
    t_next = t_n + sm * t_diff;
    
    try {
      // Try updating
      update_step(
          e_next, e_past, T_next, T_past, t_next, t_past, s_np1, s_past,
          h_np1, h_past, A_inc, E_inc, u_np1, u_past, p_np1, p_past);
    }
    // Failed adapt
    catch (const NEMLError & e) {
      nd += 1;
      if (nd >= max_divide_) {
        throw NonlinearSolverError("Exceeded maximum adaptive subdivisions");
      }
      cm /= 2;
      continue;
    }
   
    // For testing we can force subdivision
    if (force_divide_ && (nd < max_divide_)) {
      nd += 1;
      cm /= 2;
      continue;
    }

    // Tangent
    for (size_t i = 0; i < nparams()*6; i++) {
      A_old[i] += E_inc[i] * sf;
    }
    mat_mat(nparams(), 6, nparams(), A_inc, A_old, A_new);
   
    // Succeeded: advance subincrement
    cs += cm;
    std::copy(e_next, e_next+6, e_past);
    std::copy(s_np1, s_np1+6, s_past);
    std::copy(h_np1, h_np1+nstate(), h_past);
    std::copy(A_new, A_new+(nparams()*6), A_old);

    T_past = T_next;
    t_past = t_next;
    u_past = u_np1;
    p_past = p_np1;
  }

  // Extract the leading 6x6 part of the A matrix
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      A_np1[CINDEX(i,j,6)] = A_new[CINDEX(i,j,6)];
    }
  }

  delete [] h_past;
  delete [] A_inc;
  delete [] A_old;
  delete [] A_new;
  delete [] E_inc;

}

void SubstepModel_sd::update_step(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A, double * const E,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Setup the trial state
  TrialState * ts = setup(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n);

  // Take an elastic step if the model requests it
  if (elastic_step(ts, e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n)) {
    // Stress increment
    double de[6];
    sub_vec(e_np1, e_n, 6, de);

    double C[36];
    elastic_->C(T_np1, C);
    mat_vec(C, 6, de, 6, s_np1);
    for (size_t i = 0; i < 6; i++) s_np1[i] += s_n[i];
    // History
    std::copy(h_n, h_n+nstate(), h_np1);
    // Jacobian for substepping
    std::fill(A, A+(nparams()*nparams()), 0.0);
    for (size_t i = 0; i < 6; i++) A[CINDEX(i,i,nparams())] = 1.0;

    std::fill(E, E+(nparams()*6), 0.0);
    for (size_t i = 0; i < 6; i++) {
      for (size_t j = 0; j < 6; j++) {
        E[CINDEX(i,j,6)] = C[CINDEX(i,j,6)];
      }
    }

    // Energy and work
    work_and_energy(ts, e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n,
                    h_np1, h_n, u_np1, u_n, p_np1, p_n); 
    delete ts;

    return;
  }
  
  // Solve the system
  double * x = new double [nparams()]; 
  try {
    solve(this, x, ts, {rtol_, atol_, miter_, verbose_, linesearch_},
                    nullptr, A); // Keep jacobian

    // Invert the Jacobian (or idk, could go in the tangent calc)
    invert_mat(A, nparams());

    // Interpret the x vector as the updated state
    update_internal(x, e_np1, e_n, T_np1, T_n, t_np1, t_n,
                          s_np1, s_n, h_np1, h_n);

    // Get the dE matrix
    strain_partial(ts, e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n, h_np1, h_n, E);

    // Update the work and energy
    work_and_energy(ts, e_np1, e_n, T_np1, T_n, t_np1, t_n, 
                          s_np1, s_n, h_np1, h_n, u_np1, u_n,
                          p_np1, p_n);
  }
  catch (const NEMLError & e) {
    delete [] x;
    delete ts;
    throw e;
  }

  delete [] x;
  delete ts;
}

// Implementation of small strain elasticity
SmallStrainElasticity::SmallStrainElasticity(ParameterSet & params) :
    NEMLModel_sd(params)
{
  cache_history_();
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
                                          make_constant(0.0));
  pset.add_optional_parameter<bool>("truesdell", true);

  return pset;
}

std::unique_ptr<NEMLObject> SmallStrainElasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<SmallStrainElasticity>(params);
}

void SmallStrainElasticity::populate_state(History & h) const
{
}

void SmallStrainElasticity::init_state(History & h) const
{
}

void SmallStrainElasticity::update_sd_actual(
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

}

// Implementation of perfect plasticity
SmallStrainPerfectPlasticity::SmallStrainPerfectPlasticity(ParameterSet & params) :
      SubstepModel_sd(params), 
      surface_(params.get_object_parameter<YieldSurface>("surface")),
      ys_(params.get_object_parameter<Interpolate>("ys"))
{
  cache_history_();
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
                                          make_constant(0.0));
  pset.add_optional_parameter<double>("atol", 1.0e-8);
  pset.add_optional_parameter<double>("rtol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<bool>("linesearch", false);
  pset.add_optional_parameter<int>("max_divide", 4);
  pset.add_optional_parameter<bool>("force_divide", false);

  pset.add_optional_parameter<bool>("truesdell", true);

  return pset;
}

std::unique_ptr<NEMLObject> SmallStrainPerfectPlasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<SmallStrainPerfectPlasticity>(params); 
}

void SmallStrainPerfectPlasticity::populate_state(History & h) const
{

}

void SmallStrainPerfectPlasticity::init_state(History & h) const
{

}

TrialState * SmallStrainPerfectPlasticity::setup(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_n,
    const double * const h_n)
{
  SSPPTrialState * tss = new SSPPTrialState();
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, T_n,
                          s_n, h_n, *tss);
  return tss;
}

bool SmallStrainPerfectPlasticity::elastic_step(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_n,
    const double * const h_n)
{
  const SSPPTrialState * tss = static_cast<const SSPPTrialState*>(ts);
  double fv;
  surface_->f(tss->s_tr, &(tss->ys), T_np1, fv);
  return fv <= 0.0;
}

void SmallStrainPerfectPlasticity::update_internal(
    const double * const x,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n)
{
  std::copy(x, x+6, s_np1);
}

void SmallStrainPerfectPlasticity::strain_partial(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_np1, const double * const s_n,
    const double * const h_np1, const double * const h_n,
    double * const de)
{
  const SSPPTrialState * tss = static_cast<const SSPPTrialState*>(ts);
  std::fill(de, de+(6*nparams()), 0.0);

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      de[CINDEX(i,j,6)] = tss->C[CINDEX(i,j,6)];
    }
  }

}

void SmallStrainPerfectPlasticity::work_and_energy(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  const SSPPTrialState * tss = static_cast<const SSPPTrialState*>(ts);

  // Plastic work calculation
  double ds[6];
  add_vec(s_np1, s_n, 6, ds);

  double ee_np1[6];
  double S[36];
  elastic_->S(T_np1, S);
  mat_vec(S, 6, s_np1, 6, ee_np1);
  
  double dp[6];
  for (size_t i = 0; i < 6; i++) dp[i] = e_np1[i] - ee_np1[i] - tss->ep_n[i];

  p_np1 = p_n + dot_vec(ds, dp, 6) / 2.0;

  // Energy calculation (trapezoid rule)
  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  u_np1 = u_n + dot_vec(ds, de, 6) / 2.0;

}

size_t SmallStrainPerfectPlasticity::nparams() const
{
  return 7;
}

void SmallStrainPerfectPlasticity::init_x(double * const x, TrialState * ts)
{
  SSPPTrialState * tss = static_cast<SSPPTrialState *>(ts);
  std::copy(tss->s_tr, tss->s_tr+6, x);
  x[6] = 0.0;
}

void SmallStrainPerfectPlasticity::RJ(
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
  double T[6];
  for (int i=0; i<6; i++) T[i] = tss->e_np1[i] - tss->ep_n[i] - df[i] * dg;
  mat_vec(tss->C, 6, T, 6, R);
  for (int i=0; i<6; i++) R[i] = s_np1[i] - R[i];
  
  // R2
  R[6] = fv;

  // J11
  double TT[36];
  mat_mat(6, 6, 6, tss->C, ddf, TT);

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,7)] = TT[CINDEX(i,j,6)] * dg;
    }
    J[CINDEX(i,i,7)] += 1.0;
  }

  // J12
  mat_vec(tss->C, 6, df, 6, T);
  for (int i=0; i<6; i++) J[CINDEX(i,6,7)] = T[i];

  // J21
  for (int i=0; i<6; i++) J[CINDEX(6,i,7)] = df[i];

  // J22
  J[CINDEX(6,6,7)] = 0.0;
}

// Getter
double SmallStrainPerfectPlasticity::ys(double T) const 
{
  return ys_->value(T);
}

// Make this public for ease of testing
void SmallStrainPerfectPlasticity::make_trial_state(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const double * const s_n, const double * const h_n,
    SSPPTrialState & ts)
{
  ts.ys = -ys_->value(T_np1);
  ts.T = T_np1;
 
  std::copy(e_np1, e_np1+6, ts.e_np1);

  elastic_->C(T_np1, ts.C);

  double S_n[36];
  elastic_->S(T_n, S_n);
  mat_vec(S_n, 6, s_n, 6, ts.ep_n);
  for (size_t i = 0; i < 6; i++) ts.ep_n[i] = e_n[i] - ts.ep_n[i];

  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  mat_vec(ts.C, 6, de, 6, ts.s_tr);
  for (size_t i = 0; i < 6; i++) ts.s_tr[i] = s_n[i] + ts.s_tr[i];

}

// Implementation of small strain rate independent plasticity
SmallStrainRateIndependentPlasticity::SmallStrainRateIndependentPlasticity(
    ParameterSet & params) : 
      SubstepModel_sd(params),
      flow_(params.get_object_parameter<RateIndependentFlowRule>("flow"))
{
  cache_history_();
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
                                          make_constant(0.0));
  pset.add_optional_parameter<bool>("truesdell", true);

  pset.add_optional_parameter<double>("rtol", 1.0e-8);
  pset.add_optional_parameter<double>("atol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<bool>("linesearch", false);

  pset.add_optional_parameter<int>("max_divide", 4);
  pset.add_optional_parameter<bool>("force_divide", false);


  return pset;
}

std::unique_ptr<NEMLObject> SmallStrainRateIndependentPlasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<SmallStrainRateIndependentPlasticity>(params); 
}

void SmallStrainRateIndependentPlasticity::populate_state(History & hist) const
{
  flow_->set_variable_prefix(get_variable_prefix());
  flow_->populate_hist(hist);
}

void SmallStrainRateIndependentPlasticity::init_state(History & hist) const
{
  flow_->init_hist(hist);
}

TrialState * SmallStrainRateIndependentPlasticity::setup(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_n,
    const double * const h_n)
{
  SSRIPTrialState * tss = new SSRIPTrialState();
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, T_n,
                          s_n, h_n, *tss);
  return tss;
}

bool SmallStrainRateIndependentPlasticity::elastic_step(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_n,
    const double * const h_n)
{
  const SSRIPTrialState * tss = static_cast<const SSRIPTrialState *>(ts);

  double fv;
  flow_->f(tss->s_tr, &tss->h_tr[0], T_np1, fv);
  return fv <= 0.0;
}

void SmallStrainRateIndependentPlasticity::update_internal(
    const double * const x,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n)
{
  std::copy(x, x+6, s_np1);
  std::copy(x+6, x+6+nstate(), h_np1);
}

void SmallStrainRateIndependentPlasticity::strain_partial(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_np1, const double * const s_n,
    const double * const h_np1, const double * const h_n,
    double * const de)
{
  const SSRIPTrialState * tss = static_cast<const SSRIPTrialState *>(ts);
  std::fill(de, de+(6*nparams()), 0.0);

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      de[CINDEX(i,j,6)] = tss->C[CINDEX(i,j,6)];
    }
  }
}

void SmallStrainRateIndependentPlasticity::work_and_energy(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  const SSRIPTrialState * tss = static_cast<const SSRIPTrialState *>(ts);

  // Plastic work calculation
  double ds[6];
  add_vec(s_np1, s_n, 6, ds);

  double ee_np1[6];
  double S[36];
  elastic_->S(T_np1, S);
  mat_vec(S, 6, s_np1, 6, ee_np1);
  
  double  dp[6];
  for (size_t i = 0; i < 6; i++) dp[i] = e_np1[i] - ee_np1[i] - tss->ep_tr[i];

  p_np1 = p_n + dot_vec(ds, dp, 6) / 2.0;

  // Energy calculation (trapezoid rule)
  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  u_np1 = u_n + dot_vec(ds, de, 6) / 2.0;

}

size_t SmallStrainRateIndependentPlasticity::nparams() const
{
  // My convention: plastic strain, history, consistency parameter
  return 6 + nstate() + 1;
}

void SmallStrainRateIndependentPlasticity::init_x(double * const x, TrialState * ts)
{
  SSRIPTrialState * tss = static_cast<SSRIPTrialState *>(ts);
  std::copy(tss->s_tr, tss->s_tr+6, x);
  std::copy(&tss->h_tr[0], &tss->h_tr[0]+nstate(), &x[6]);
  x[6+nstate()] = 0.0; // consistency parameter
}

void SmallStrainRateIndependentPlasticity::RJ(const double * const x, 
                                             TrialState * ts, 
                                             double * const R, double * const J)
{
  SSRIPTrialState * tss = static_cast<SSRIPTrialState *>(ts);

  // Again no idea why the compiler can't do this
  int nh = nstate();

  // Setup from current state
  const double * const s_np1 = &x[0];
  const double * const alpha  = &x[6];
  const double & dg = x[6+nh];

  // Residual calculation
  double g[6];
  flow_->g(s_np1, alpha, tss->T, g); 
  std::vector<double> hv(nh);
  double * h = &hv[0];
  flow_->h(s_np1, alpha, tss->T, h);
  double f;
  flow_->f(s_np1, alpha, tss->T, f);

  double R1[6];
  for (int i=0; i<6; i++) {
    R1[i] = tss->e_np1[i] - tss->ep_tr[i] - g[i] * dg;
  }
  mat_vec(tss->C, 6, R1, 6, R);
  for (int i=0; i<6; i++) {
    R[i] = s_np1[i] - R[i];
  }

  for (int i=0; i<nh; i++) {
    R[i+6] = alpha[i] - tss->h_tr[i] - h[i] * dg;
  }
  R[6+nh] = f;

  // Now the jacobian calculation...
  int n = nparams();
  
  // J11
  double gs[36];
  double J11[36];
  flow_->dg_ds(s_np1, alpha, tss->T, gs);
  mat_mat(6, 6, 6, tss->C, gs, J11);
  for (int i=0; i<36; i++) J11[i] = J11[i] * dg;
  for (int i=0; i<6; i++) J11[CINDEX(i,i,6)] += 1.0;

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,n)] = J11[CINDEX(i,j,6)];
    }
  }
  
  // J12
  std::vector<double> gav(6*nh);
  double * ga = &gav[0];
  std::vector<double> J12v(6*nh);
  double * J12 = &J12v[0];
  flow_->dg_da(s_np1, alpha, tss->T, ga);
  mat_mat(6, nh, 6, tss->C, ga, J12);

  for (int i=0; i<6*nh; i++) J12[i] *= dg;
  for (int i=0; i<6; i++) {
    for (int j=0; j < nh; j++) {
      J[CINDEX(i,(j+6),n)] = J12[CINDEX(i,j,nh)];
    }
  }

  // J13
  double J13[6];
  mat_vec(tss->C, 6, g, 6, J13);
  for (int i=0; i<6; i++) {
    J[CINDEX(i,(6+nh),n)] = J13[i];
  }

  // J21
  std::vector<double> J21v(nh*6);
  double * J21 = &J21v[0];
  flow_->dh_ds(s_np1, alpha, tss->T, J21);
  for (int i=0; i<nh*6; i++) J21[i] = J21[i] * dg;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,n)] = -J21[CINDEX(i,j,6)];
    }
  }

  // J22
  std::vector<double> J22v(nh*nh);
  double * J22 = &J22v[0];
  flow_->dh_da(s_np1, alpha, tss->T, J22);
  for (int i=0; i<nh*nh; i++) J22[i] *= dg;
  for (int i=0; i<nh; i++) J22[CINDEX(i,i,nh)] -= 1.0;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<nh; j++) {
      J[CINDEX((i+6), (j+6), n)] = -J22[CINDEX(i,j,nh)];
    }
  }

  // J23
  for (int i=0; i<nh; i++) {
    J[CINDEX((i+6), (6+nh), n)] = -h[i];
  }

  // J31
  double J31[6];
  flow_->df_ds(s_np1, alpha, tss->T, J31);
  for (int i=0; i<6; i++) {
    J[CINDEX((6+nh), i, n)] = J31[i];
  }

  // J32
  std::vector<double> J32v(nh);
  double * J32 = &J32v[0];
  flow_->df_da(s_np1, alpha, tss->T, J32);
  for (int i=0; i<nh; i++) {
    J[CINDEX((6+nh), (i+6), n)] = J32[i];
  }

  // J33
  J[CINDEX((6+nh), (6+nh), n)] = 0.0;
}

const std::shared_ptr<const LinearElasticModel> SmallStrainRateIndependentPlasticity::elastic() const
{
  return elastic_;
}

void SmallStrainRateIndependentPlasticity::make_trial_state(
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

  ts.h_tr.resize(nstate());
  std::copy(h_n, h_n+nstate(), ts.h_tr.begin());
  // Calculate the trial stress
  double ee[6];
  sub_vec(e_np1, ts.ep_tr, 6, ee);
  elastic_->C(T_np1, ts.C);

  mat_vec(ts.C, 6, ee, 6, ts.s_tr);
  // Store temp
  ts.T = T_np1;
}

// Implement creep + plasticity
// Implementation of small strain rate independent plasticity
//
SmallStrainCreepPlasticity::SmallStrainCreepPlasticity(
    ParameterSet & params) :
      NEMLModel_sd(params),
      plastic_(params.get_object_parameter<NEMLModel_sd>("plastic")),
      creep_(params.get_object_parameter<CreepModel>("creep")),
      rtol_(params.get_parameter<double>("rtol")),
      atol_(params.get_parameter<double>("atol")),
      sf_(params.get_parameter<double>("sf")),
      miter_(params.get_parameter<int>("miter")),
      verbose_(params.get_parameter<bool>("verbose")),
      linesearch_(params.get_parameter<bool>("linesearch"))
{
  cache_history_();
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
                                           make_constant(0.0));
  pset.add_optional_parameter<double>("atol", 1.0e-10);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<double>("sf", 1.0e6);
  pset.add_optional_parameter<double>("rtol", 1.0e-12);
  pset.add_optional_parameter<bool>("linesearch", true); // Note this

  pset.add_optional_parameter<bool>("truesdell", true);

  return pset;
}

std::unique_ptr<NEMLObject> SmallStrainCreepPlasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<SmallStrainCreepPlasticity>(params); 
}

void SmallStrainCreepPlasticity::populate_state(History & hist) const
{
  plastic_->set_variable_prefix(get_variable_prefix());
  plastic_->populate_state(hist);
  hist.add<Symmetric>(prefix("plastic_strain"));
}

void SmallStrainCreepPlasticity::init_state(History & hist) const
{
  hist.get<Symmetric>(prefix("plastic_strain")) = Symmetric::zero();
  plastic_->init_state(hist);
}

void SmallStrainCreepPlasticity::update_sd_actual(
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

  std::vector<double> xv(nparams());
  double * x = &xv[0];
  solve(this, x, &ts, {rtol_, atol_, miter_, verbose_, linesearch_});

  // Store the ep strain
  std::copy(x, x+6, h_np1);

  // Do the plastic update to get the new history and stress
  double A[36];
  plastic_->update_sd_actual(x, ts.ep_strain, T_np1, T_n,
                             t_np1, t_n, s_np1, s_n,
                             &h_np1[6], &h_n[6],
                             A, u_np1, u_n, p_np1, p_n);

  // Do the creep update to get a tangent component
  double creep_old[6];
  double creep_new[6];
  double B[36];
  for (int i=0; i<6; i++) {
    creep_old[i] = e_n[i] - ts.ep_strain[i];
  }
  creep_->update(s_np1, creep_new, creep_old, T_np1, T_n,
                 t_np1, t_n, B);

  // Form the relatively simple tangent
  form_tangent_(A, B, A_np1);

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
  
}

size_t SmallStrainCreepPlasticity::nparams() const
{
  // Just the elastic-plastic strain
  return 6;
}

void SmallStrainCreepPlasticity::init_x(double * const x, TrialState * ts)
{
  SSCPTrialState * tss = static_cast<SSCPTrialState*>(ts);

  // Start out at last step's value
  std::copy(tss->ep_strain, tss->ep_strain + 6, x);
}

void SmallStrainCreepPlasticity::RJ(const double * const x, TrialState * ts, 
                                   double * const R, double * const J)
{
  SSCPTrialState * tss = static_cast<SSCPTrialState*>(ts);

  // First update the elastic-plastic model
  double s_np1[6];
  double A_np1[36];
  std::vector<double> h_np1;
  h_np1.resize(plastic_->nstate());
  double u_np1, u_n;
  double p_np1, p_n;
  u_n = 0.0;
  p_n = 0.0;

  double * hist = (h_np1.empty() ? nullptr : &h_np1[0]);
  double * hist_tss = (tss->h_n.empty() ? nullptr : &(tss->h_n[0]));

  plastic_->update_sd_actual(x, tss->ep_strain, tss->T_np1, tss->T_n,
                      tss->t_np1, tss->t_n, s_np1, tss->s_n,
                      hist, hist_tss, A_np1,
                      u_np1, u_n, p_np1, p_n);

  // Then update the creep strain
  double creep_old[6];
  double creep_new[6];
  double B[36];
  for (int i=0; i<6; i++) {
    creep_old[i] = tss->e_n[i] - tss->ep_strain[i];
  }
  creep_->update(s_np1, creep_new, creep_old, tss->T_np1, tss->T_n,
                 tss->t_np1, tss->t_n, B);

  // Form the residual
  for (int i=0; i<6; i++) {
    R[i] = (x[i] + creep_new[i] - tss->e_np1[i]) * sf_;
  }
  
  // The Jacobian is a straightforward combination of the two derivatives
  mat_mat(6, 6, 6, B, A_np1, J);
  for (int i=0; i<6; i++) J[CINDEX(i,i,6)] += 1.0;
  for (int i=0; i<36; i++) J[i] *= sf_;
}

void SmallStrainCreepPlasticity::make_trial_state(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const double * const s_n, const double * const h_n,
    SSCPTrialState & ts)
{
  int nh = plastic_->nstate();
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

}

void SmallStrainCreepPlasticity::form_tangent_(
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

}

void SmallStrainCreepPlasticity::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  plastic_->set_elastic_model(emodel);
}

// Start general integrator implementation
GeneralIntegrator::GeneralIntegrator(ParameterSet & params) :
    SubstepModel_sd(params),
    rule_(params.get_object_parameter<GeneralFlowRule>("rule")),
    skip_first_(params.get_parameter<bool>("skip_first_step"))
{
  cache_history_();
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
                                          make_constant(0.0));
  pset.add_optional_parameter<bool>("truesdell", true);
  
  pset.add_optional_parameter<double>("rtol", 1.0e-8);
  pset.add_optional_parameter<double>("atol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<bool>("linesearch", false);
  pset.add_optional_parameter<int>("max_divide", 4);
  pset.add_optional_parameter<bool>("force_divide", false);
  pset.add_optional_parameter<bool>("skip_first_step", false);

  return pset;
}

std::unique_ptr<NEMLObject> GeneralIntegrator::initialize(ParameterSet & params)
{
  return neml::make_unique<GeneralIntegrator>(params); 
}

TrialState * GeneralIntegrator::setup(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_n,
    const double * const h_n)
{
  GITrialState * tss = new GITrialState();
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, 
                   s_n, h_n, *tss);
  return tss;
}

bool GeneralIntegrator::elastic_step(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_n,
    const double * const h_n)
{
  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  
  double dt = t_np1 - t_n;

  if (dt < std::numeric_limits<double>::epsilon()) {
    return true;
  }

  return false;
}

void GeneralIntegrator::update_internal(
    const double * const x,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n)
{
  std::copy(x, x+6, s_np1);
  std::copy(x+6, x+6+nstate(), h_np1);

}

void GeneralIntegrator::strain_partial(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const double * const s_np1, const double * const s_n,
    const double * const h_np1, const double * const h_n,
    double * de)
{
  const GITrialState * tss = static_cast<const GITrialState*>(ts);
  
  double * estress = new double [6*6];

  rule_->ds_de(s_np1, h_np1, tss->e_dot, tss->T, tss->Tdot, estress);
  
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      de[CINDEX(i,j,6)] = estress[CINDEX(i,j,6)];
    }
  }

  delete [] estress;

  double * ehist = new double [6*nstate()];

  rule_->da_de(s_np1, h_np1, tss->e_dot, tss->T, tss->Tdot, ehist);
  for (size_t i = 0; i < nstate(); i++) {
    for (size_t j = 0; j < 6; j++) {
      de[CINDEX((i+6),j,6)] = ehist[CINDEX(i,j,6)];
    }
  }

  delete [] ehist;

}

void GeneralIntegrator::work_and_energy(
    const TrialState * ts,
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  const GITrialState * tss = static_cast<const GITrialState *>(ts);

  // Energy calculation (trapezoid rule)
  double de[6];
  double ds[6];
  sub_vec(e_np1, e_n, 6, de);
  add_vec(s_np1, s_n, 6, ds);
  for (int i=0; i<6; i++) ds[i] /= 2.0;
  u_np1 = u_n + dot_vec(ds, de, 6);

  // Need a special call
  double p_dot_np1;
  rule_->work_rate(s_np1, h_np1, tss->e_dot, T_np1, tss->Tdot, p_dot_np1);
  double p_dot_n;
  rule_->work_rate(s_n, h_n, tss->e_dot, T_n, tss->Tdot, p_dot_n);
  p_np1 = p_n + (p_dot_np1 + p_dot_n)/2.0 * tss->dt;

}

void GeneralIntegrator::populate_state(History & hist) const
{
  rule_->set_variable_prefix(get_variable_prefix());
  rule_->populate_hist(hist);
}

void GeneralIntegrator::init_state(History & hist) const
{
  rule_->init_hist(hist);
}

size_t GeneralIntegrator::nparams() const
{
  return 6 + nstate();
}

void GeneralIntegrator::init_x(double * const x, TrialState * ts)
{
  GITrialState * tss = static_cast<GITrialState*>(ts);
  std::copy(tss->s_guess, tss->s_guess+6, x);
  std::copy(tss->h_n.begin(), tss->h_n.end(), &x[6]);

  rule_->override_guess(x);
}

void GeneralIntegrator::RJ(const double * const x, TrialState * ts,
                          double * const R, double * const J)
{
  GITrialState * tss = static_cast<GITrialState*>(ts);

  // Setup
  const double * s_np1 = x;
  const double * const h_np1 = &x[6];
  
  // Helps with vectorization
  // Really as I declared both const this shouldn't be necessary but hey
  // I don't design optimizing compilers for a living
  int nstate = this->nstate();
  int nparams = this->nparams();

  // Residual calculation
  rule_->s(s_np1, h_np1, tss->e_dot, tss->T, tss->Tdot, R);
  for (int i=0; i<6; i++) {
    R[i] = s_np1[i] - tss->s_n[i] - R[i] * tss->dt;
  }
  rule_->a(s_np1, h_np1, tss->e_dot, tss->T, tss->Tdot, &R[6]);
  for (int i=0; i<nstate; i++) {
    R[i+6] = h_np1[i] - tss->h_n[i] - R[i+6] * tss->dt;
  }

  // Jacobian calculation
  double J11[36];
  rule_->ds_ds(s_np1, h_np1, tss->e_dot, tss->T, tss->Tdot, J11);
  for (int i=0; i<36; i++) J11[i] *= tss->dt;
  for (int i=0; i<6; i++) J11[CINDEX(i,i,6)] -= 1.0;
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,nparams)] = -J11[CINDEX(i,j,6)];
    }
  }
  
  std::vector<double> J12v(6*nstate);
  double * J12 = &J12v[0];
  rule_->ds_da(s_np1, h_np1, tss->e_dot, tss->T, tss->Tdot, J12);
  for (int i=0; i<6; i++) {
    for (int j=0; j<nstate; j++) {
      J[CINDEX(i,(j+6),nparams)] = -J12[CINDEX(i,j,nstate)] * tss->dt;
    }
  }
  
  std::vector<double> J21v(nstate*6);
  double * J21 = &J21v[0];
  rule_->da_ds(s_np1, h_np1, tss->e_dot, tss->T, tss->Tdot, J21);
  for (int i=0; i<nstate; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,nparams)] = -J21[CINDEX(i,j,6)] * tss->dt;
    }
  }
  
  std::vector<double> J22v(nstate*nstate);
  double * J22 = &J22v[0];
  rule_->da_da(s_np1, h_np1, tss->e_dot, tss->T, tss->Tdot, J22);

  // More vectorization
  double dt = tss->dt;
  for (int i=0; i<nstate*nstate; i++) J22[i] *= dt;
  for (int i=0; i<nstate; i++) J22[CINDEX(i,i,nstate)] -= 1.0;

  for (int i=0; i<nstate; i++) {
    for (int j=0; j<nstate; j++) {
      J[CINDEX((i+6),(j+6),nparams)] = -J22[CINDEX(i,j,nstate)];
    }
  }
}


void GeneralIntegrator::make_trial_state(
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
  
  // Last stress
  std::copy(s_n, s_n+6, ts.s_n);

  // Last history
  ts.h_n.resize(nstate());
  std::copy(h_n, h_n+nstate(), ts.h_n.begin());

  // Elastic guess
  double C[36];
  elastic_->C(T_np1, C);
  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  mat_vec(C, 6, de, 6, ts.s_guess);
  add_vec(ts.s_guess, s_n, 6, ts.s_guess);

  // Special logic...
  if ((t_n == 0.0) && (skip_first_))
    std::copy(s_n, s_n+6, ts.s_guess);

}

void GeneralIntegrator::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  rule_->set_elastic_model(emodel);

}

// Start KMRegimeModel
KMRegimeModel::KMRegimeModel(ParameterSet & params) :
    NEMLModel_sd(params), 
    models_(params.get_object_parameter_vector<NEMLModel_sd>("models")), 
    gs_(params.get_parameter<std::vector<double>>("gs")),
    kboltz_(params.get_parameter<double>("kboltz")), 
    b_(params.get_parameter<double>("b")), 
    eps0_(params.get_parameter<double>("eps0"))
{
  cache_history_();
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
                                          make_constant(0.0));

  pset.add_optional_parameter<bool>("truesdell", true);

  return pset;
}

std::unique_ptr<NEMLObject> KMRegimeModel::initialize(ParameterSet & params)
{
  return neml::make_unique<KMRegimeModel>(params); 
}

void KMRegimeModel::update_sd_actual(
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
      return models_[i]->update_sd_actual(e_np1, e_n, T_np1, T_n, t_np1, t_n, 
                                   s_np1, s_n, h_np1, h_n, A_np1, u_np1, u_n,
                                   p_np1, p_n);
    }
  }
  return models_.back()->update_sd_actual(e_np1, e_n, T_np1, T_n, t_np1, t_n,
                                   s_np1, s_n, h_np1, h_n, A_np1, u_np1, u_n,
                                   p_np1, p_n);
}

void KMRegimeModel::populate_state(History & hist) const
{
  for (auto model : models_)
    model->set_variable_prefix(get_variable_prefix());

  return models_[0]->populate_hist(hist);
}

void KMRegimeModel::init_state(History & hist) const
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

void KMRegimeModel::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    (*it)->set_elastic_model(emodel);
  }
}

} // namespace neml
