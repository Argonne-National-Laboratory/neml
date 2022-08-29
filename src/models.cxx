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

void NEMLModel::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Wrap everything and just call the interface
  Symmetric E_np1(e_np1);
  Symmetric E_n(e_n);

  Symmetric S_np1(s_np1);
  Symmetric S_n(s_n);

  History H_np1 = gather_history_(h_np1);
  History H_n = gather_history_(h_n);

  SymSymR4 A(A_np1);

  update_sd_interface(E_np1, E_n, T_np1, T_n, t_np1, t_n, S_np1, S_n,
                      H_np1, H_n, A, u_np1, u_n, p_np1, p_n);
}

void NEMLModel::update_ld_inc(
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
  // Wrap everything and just call the interface
  Symmetric D_np1(d_np1);
  Symmetric D_n(d_n);
  Skew W_np1(w_np1);
  Skew W_n(w_n);
  Symmetric S_np1(s_np1);
  Symmetric S_n(s_n);
  SymSymR4 A(A_np1);
  SymSkewR4 B(B_np1);

  History H_np1 = gather_history_(h_np1);
  History H_n = gather_history_(h_n);

  update_ld_inc_interface(D_np1, D_n, W_np1, W_n, T_np1, T_n, t_np1, t_n, 
                S_np1, S_n, H_np1, H_n, A, B, u_np1, u_n, 
                p_np1, p_n);
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

void NEMLModel::elastic_strains(const double * const s_np1,
                                double T_np1, const double * const h_np1,
                                double * const e_np1) const
{
  Symmetric S_np1(s_np1);
  History H_np1 = gather_history_(h_np1);
  Symmetric E_np1(e_np1);
  E_np1 = elastic_strains_interface(S_np1, T_np1, H_np1);
}

Symmetric NEMLModel::elastic_strains_interface(const Symmetric & s_np1, double T_np1, const History & h_np1) const
{
  return Symmetric();
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
  populate_hist(stored_hist_);
  populate_state(stored_state_);
  populate_static(stored_static_);
}

History NEMLModel::gather_history_(double * data) const
{
  History h = gather_blank_history_();
  h.set_data(data);
  return h;
}

History NEMLModel::gather_history_(const double * data) const
{
  History h = gather_blank_history_();
  h.set_data(const_cast<double*>(data));
  return h;
}

History NEMLModel::gather_blank_history_() const
{
  return stored_hist_;
}

History NEMLModel::gather_state_(double * data) const
{
  History h = gather_blank_state_();
  h.set_data(data);
  return h;
}

History NEMLModel::gather_state_(const double * data) const
{
  History h = gather_blank_state_();
  h.set_data(const_cast<double*>(data));
  return h;
}

History NEMLModel::gather_blank_state_() const
{
  return stored_state_;
}

// NEMLModel_sd implementation
NEMLModel_sd::NEMLModel_sd(ParameterSet & params) :
      NEMLModel(params), 
      elastic_(params.get_object_parameter<LinearElasticModel>("elastic")), 
      alpha_(params.get_object_parameter<Interpolate>("alpha")),
      truesdell_(params.get_parameter<bool>("truesdell"))
{

}

void NEMLModel_sd::update_ld_inc_interface(
    const Symmetric & d_np1, const Symmetric & d_n, 
    const Skew & w_np1, const Skew & w_n, 
    double T_np1, double T_n, 
    double t_np1, double t_n,
    Symmetric & s_np1, const Symmetric & s_n,
    History & h_np1, const History & h_n,
    SymSymR4 & A_np1, SymSkewR4 & B_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  Symmetric D = d_np1 - d_n;
  Skew W = w_np1 - w_n;
  if (not truesdell_)
    D = Symmetric::zero();
  
  SymSymR4 base_A_np1;
  Symmetric stress_prime;
  Symmetric stress_prime_n = h_n.get<Symmetric>(prefix("small_stress"));

  update_sd_interface(
      d_np1, d_n, 
      T_np1, T_n, t_np1, t_n, 
      stress_prime,
      stress_prime_n,
      h_np1, h_n,
      base_A_np1, u_np1, u_n, p_np1, p_n);
  
  h_np1.get<Symmetric>(prefix("small_stress")) = stress_prime;
  Symmetric dS = stress_prime - 
      h_n.get<Symmetric>(prefix("small_stress"));

  s_np1 = truesdell_update_sym(D, W, s_n, dS);
  calc_tangent_(D, W, base_A_np1, s_np1, A_np1, B_np1);
}

void NEMLModel_sd::update_sd_interface(
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    Symmetric & s_np1, Symmetric & s_n,
    History & h_np1, const History & h_n,
    SymSymR4 & A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  auto [H_np1, F_np1] = split_state(h_np1);
  auto [H_n, F_n] = split_state(h_n);
  
  update_sd_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n,
            H_np1, H_n, 
            A_np1, u_np1, u_n, p_np1, p_n);
}

void NEMLModel_sd::init_static(History & h) const
{
  h.get<Symmetric>(prefix("small_stress")) = Symmetric::zero();
}

void NEMLModel_sd::populate_static(History & h) const
{
  // This is a legacy from the old flat vector system.  Composite models
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

Symmetric NEMLModel_sd::elastic_strains_interface(const Symmetric & s_np1, double T_np1, const History & h_np1) const
{
  return elastic_->S(T_np1).dot(s_np1);
}

void NEMLModel_sd::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
}

void NEMLModel_sd::calc_tangent_(const Symmetric & D, const Skew & W,
                                 const SymSymR4 & C,  const Symmetric & S,
                                 SymSymR4 & A, SymSkewR4 & B)
{
  double J[81];
  double O[81];

  truesdell_mat(D.data(), W.data(), J);
  invert_mat(J, 9);

  truesdell_tangent_outer(S.data(), O);

  double F[81];
  mandel2full(C.data(), F);
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

  full2mandel(Dpart, A.s());
  full2skew(Wpart, B.s());

  // IMPORTANT TODO: go back and find where you dropped the factor of 2...
  for (int i=0; i<18; i++) B.s()[i] *= 2.0;

}

// NEMLModel_ldi implementation
NEMLModel_ldi::NEMLModel_ldi(ParameterSet & params) :
      NEMLModel(params)
{

}

void NEMLModel_ldi::update_sd_interface(
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    Symmetric & s_np1, Symmetric & s_n,
    History & h_np1, const History & h_n,
    SymSymR4 & A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  Skew W;
  SymSkewR4 B;
  update_ld_inc_interface(e_np1, e_n, W, W, T_np1, T_n, t_np1, t_n,
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

void SubstepModel_sd::update_sd_state(
    const Symmetric & E_np1, const Symmetric & E_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    Symmetric & S_np1, const Symmetric & S_n,
    History & H_np1, const History & H_n,
    SymSymR4 & AA_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Setup the substep parameters
  int nd = 0;                     // Number of times we subdivided
  int tf = pow(2, max_divide_);   // Total integer step count
  int cm = tf;                    // Attempted integer step (sf = cm/tf)
  int cs = 0;                     // Accumulated steps
  
  // Strain, time, and temperature increments
  Symmetric E_diff = E_np1 - E_n;
  double T_diff = T_np1 - T_n;
  double t_diff = t_np1 - t_n;

  // Previous subincrement quantities
  Symmetric E_past = E_n;
  Symmetric S_past = S_n;
  History H_past = H_n.deepcopy(); // We're not allowed to wrote to H_n, so need copy
  double T_past = T_n;
  double t_past = t_n;
  double u_past = u_n;
  double p_past = p_n;

  // Targets
  Symmetric E_next;
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
    E_next = E_n + sm * E_diff;
    T_next = T_n + sm * T_diff;
    t_next = t_n + sm * t_diff;
    
    try {
      // Try updating
      update_step(
          E_next, E_past, T_next, T_past, t_next, t_past, S_np1, S_past,
          H_np1, H_past, A_inc, E_inc, u_np1, u_past, p_np1, p_past);
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
    E_past = E_next;
    S_past = S_np1;
    H_past = H_np1;
    std::copy(A_new, A_new+(nparams()*6), A_old);

    T_past = T_next;
    t_past = t_next;
    u_past = u_np1;
    p_past = p_np1;
  }

  // Extract the leading 6x6 part of the A matrix
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      AA_np1(i,j) = A_new[CINDEX(i,j,6)];
    }
  }

  delete [] A_inc;
  delete [] A_old;
  delete [] A_new;
  delete [] E_inc;
}

void SubstepModel_sd::update_step(
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    Symmetric & s_np1, const Symmetric & s_n,
    History & h_np1, const History & h_n,
    double * const A, double * const E,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Setup the trial state
  std::unique_ptr<TrialState> ts = setup(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n);

  // Take an elastic step if the model requests it
  if (elastic_step(ts.get(), e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n)) {
    // Strain increment
    Symmetric de = e_np1 - e_n;

    // Stress increment
    SymSymR4 C = elastic_->C(T_np1);
    s_np1 = C.dot(de) + s_n;
    
    // History
    h_np1 = h_n;

    // Jacobian for substepping
    std::fill(A, A+(nparams()*nparams()), 0.0);
    for (size_t i = 0; i < 6; i++) A[CINDEX(i,i,nparams())] = 1.0;

    std::fill(E, E+(nparams()*6), 0.0);
    for (size_t i = 0; i < 6; i++) {
      for (size_t j = 0; j < 6; j++) {
        E[CINDEX(i,j,6)] = C(i,j);
      }
    }

    // Energy and work
    work_and_energy(ts.get(), e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n,
                    h_np1, h_n, u_np1, u_n, p_np1, p_n); 
    return;
  }
  
  // Solve the system
  std::vector<double> x(nparams());
  try {
    solve(this, &x[0], ts.get(), {rtol_, atol_, miter_, verbose_, linesearch_},
                    nullptr, A); // Keep jacobian

    // Invert the Jacobian (or idk, could go in the tangent calc)
    invert_mat(A, nparams());

    // Interpret the x vector as the updated state
    update_internal(&x[0], e_np1, e_n, T_np1, T_n, t_np1, t_n,
                          s_np1, s_n, h_np1, h_n);

    // Get the dE matrix
    strain_partial(ts.get(), e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n, 
                   h_np1, h_n, E);

    // Update the work and energy
    work_and_energy(ts.get(), e_np1, e_n, T_np1, T_n, t_np1, t_n, 
                          s_np1, s_n, h_np1, h_n, u_np1, u_n,
                          p_np1, p_n);
  }
  catch (const NEMLError & e) {
    throw e;
  }
}

void SubstepModel_sd::work_and_energy(
    const TrialState * ts,
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const Symmetric & s_np1, const Symmetric & s_n,
    const History & h_np1, const History & h_n,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  Symmetric ep_np1 = e_np1 - elastic_->S(T_np1).dot(s_np1);
  Symmetric ep_n = e_n - elastic_->S(T_n).dot(s_n);

  auto [du, dp] = trapezoid_energy(e_np1, e_n, ep_np1, ep_n,
                                   s_np1, s_n);

  u_np1 = u_n + du;
  p_np1 = p_n + dp;
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

void SmallStrainElasticity::update_sd_state(
    const Symmetric & E_np1, const Symmetric & E_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    Symmetric & S_np1, const Symmetric & S_n,
    History & H_np1, const History & H_n,
    SymSymR4 & AA_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  AA_np1 = elastic_->C(T_np1);
  S_np1 = AA_np1.dot(E_np1);

  auto [du, dp] = trapezoid_energy(E_np1, E_n,
                                   Symmetric::zero(),
                                   Symmetric::zero(),
                                   S_np1, S_n);
  u_np1 = u_n + du;
  p_np1 = p_n + dp;
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

std::unique_ptr<TrialState> SmallStrainPerfectPlasticity::setup(
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const Symmetric & s_n,
    const History & h_n)
{
  Symmetric ep_n = e_n - elastic_->S(T_n).dot(s_n);
  Symmetric s_tr = s_n + elastic_->C(T_np1).dot(e_np1 - e_n);
  return std::make_unique<SSPPTrialState>(
      e_np1, ep_n, s_tr, elastic_->C(T_np1), -ys_->value(T_np1),
      T_np1);
}

bool SmallStrainPerfectPlasticity::elastic_step(
    const TrialState * ts, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n, const Symmetric & s_n,
    const History & h_n)
{
  const SSPPTrialState * tss = static_cast<const SSPPTrialState*>(ts);
  double fv;
  surface_->f(tss->s_tr.data(), &(tss->ys), T_np1, fv);
  return fv <= 0.0;
}

void SmallStrainPerfectPlasticity::update_internal(
    const double * const x, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    Symmetric & s_np1, const Symmetric & s_n,
    History & h_np1, const History & h_n)
{
  s_np1.copy_data(x);
}

void SmallStrainPerfectPlasticity::strain_partial(
    const TrialState * ts, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const Symmetric & s_np1, const Symmetric & s_n, 
    const History & h_np1, const History & h_n, double * de)
{
  const SSPPTrialState * tss = static_cast<const SSPPTrialState*>(ts);
  std::fill(de, de+(6*nparams()), 0.0);

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      de[CINDEX(i,j,6)] = tss->C(i,j);
    }
  }
}

size_t SmallStrainPerfectPlasticity::nparams() const
{
  return 7;
}

void SmallStrainPerfectPlasticity::init_x(double * const x, TrialState * ts)
{
  SSPPTrialState * tss = static_cast<SSPPTrialState *>(ts);
  std::copy(tss->s_tr.data(), tss->s_tr.data()+6, x);
  x[6] = 0.0;
}

void SmallStrainPerfectPlasticity::RJ(
    const double * const x, TrialState * ts, double * const R,
    double * const J)
{
  SSPPTrialState * tss = static_cast<SSPPTrialState *>(ts);
  Symmetric s_np1(x);
  double dg = x[6];

  // Common things
  double fv;
  surface_->f(s_np1.data(), &tss->ys, tss->T, fv);
  
  Symmetric df;
  surface_->df_ds(s_np1.data(), &tss->ys, tss->T, df.s());

  SymSymR4 ddf;
  surface_->df_dsds(s_np1.data(), &tss->ys, tss->T, ddf.s());

  // R1
  Symmetric R1 = s_np1 - tss->C.dot(tss->e_np1 - tss->ep_n - df * dg);
  for (size_t i = 0; i < 6; i++) R[i] = R1(i);
  
  // R2
  R[6] = fv;

  // J11
  SymSymR4 J11 = tss->C.dot(ddf)*dg + SymSymR4::id();

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,7)] = J11(i,j);
    }
  }

  // J12
  Symmetric J12 = tss->C.dot(df);
  for (int i=0; i<6; i++) J[CINDEX(i,6,7)] = J12(i);

  // J21
  for (int i=0; i<6; i++) J[CINDEX(6,i,7)] = df(i);

  // J22
  J[CINDEX(6,6,7)] = 0.0;
}

// Getter
double SmallStrainPerfectPlasticity::ys(double T) const 
{
  return ys_->value(T);
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

std::unique_ptr<TrialState> SmallStrainRateIndependentPlasticity::setup(
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const Symmetric & s_n,
    const History & h_n)
{
  Symmetric ep_tr = e_n - elastic_->S(T_n).dot(s_n);
  SymSymR4 C = elastic_->C(T_np1);
  Symmetric s_tr = C.dot(e_np1 - ep_tr);

  return std::make_unique<SSRIPTrialState>(e_np1, ep_tr, s_tr, C,
                                           h_n, T_np1);
}

bool SmallStrainRateIndependentPlasticity::elastic_step(
    const TrialState * ts, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n, const Symmetric & s_n,
    const History & h_n)
{
  const SSRIPTrialState * tss = static_cast<const SSRIPTrialState *>(ts);

  double fv;
  flow_->f(tss->s_tr.data(), tss->h_tr.rawptr(), T_np1, fv);
  return fv <= 0.0;
}

void SmallStrainRateIndependentPlasticity::update_internal(
    const double * const x, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    Symmetric & s_np1, const Symmetric & s_n,
    History & h_np1, const History & h_n)
{
  s_np1.copy_data(x);
  h_np1.copy_data(x+6);
}

void SmallStrainRateIndependentPlasticity::strain_partial(
    const TrialState * ts, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const Symmetric & s_np1, const Symmetric & s_n, 
    const History & h_np1, const History & h_n, double * de)
{
  const SSRIPTrialState * tss = static_cast<const SSRIPTrialState *>(ts);
  std::fill(de, de+(6*nparams()), 0.0);

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      de[CINDEX(i,j,6)] = tss->C(i,j);
    }
  }
}

size_t SmallStrainRateIndependentPlasticity::nparams() const
{
  // My convention: plastic strain, history, consistency parameter
  return 6 + nstate() + 1;
}

void SmallStrainRateIndependentPlasticity::init_x(double * const x, TrialState * ts)
{
  SSRIPTrialState * tss = static_cast<SSRIPTrialState *>(ts);
  std::copy(tss->s_tr.data(), tss->s_tr.data()+6, x);
  std::copy(tss->h_tr.rawptr(), tss->h_tr.rawptr()+nstate(), &x[6]);
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
  Symmetric s_np1(x);
  History alpha = gather_state_(&x[6]);
  double dg = x[6+nh];

  // Residual calculation
  Symmetric g;
  flow_->g(s_np1.data(), alpha.rawptr(), tss->T, g.s()); 
  History h = gather_blank_state_();
  flow_->h(s_np1.data(), alpha.rawptr(), tss->T, h.rawptr());
  double f;
  flow_->f(s_np1.data(), alpha.rawptr(), tss->T, f);
  
  // R1
  Symmetric R1 = s_np1 - tss->C.dot(tss->e_np1 - tss->ep_tr - g * dg);
  for (size_t i = 0; i<6; i++) R[i] = R1(i);

  // R2
  for (int i=0; i<nh; i++) {
    R[i+6] = alpha.rawptr()[i] - tss->h_tr.rawptr()[i] - h.rawptr()[i] * dg;
  }
  R[6+nh] = f;

  // Now the jacobian calculation...
  int n = nparams();
  
  // J11
  SymSymR4 gs, J11;
  flow_->dg_ds(s_np1.data(), alpha.rawptr(), tss->T, gs.s());
  J11 = tss->C.dot(gs) * dg + SymSymR4::id();
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,n)] = J11(i,j);
    }
  }
  
  // J12 (transpose?)
  History ga = gather_blank_state_().derivative<Symmetric>();
  flow_->dg_da(s_np1.data(), alpha.rawptr(), tss->T, ga.rawptr());
  History J12 = ga.premultiply(tss->C) * dg;
  for (int i=0; i<6; i++) {
    for (int j=0; j < nh; j++) {
      J[CINDEX(i,(j+6),n)] = J12.rawptr()[CINDEX(i,j,nh)];
    }
  }

  // J13
  Symmetric J13 = tss->C.dot(g);
  for (int i=0; i<6; i++) {
    J[CINDEX(i,(6+nh),n)] = J13.data()[i];
  }

  // J21 (transpose?)
  History J21 = gather_blank_state_().derivative<Symmetric>();
  flow_->dh_ds(s_np1.data(), alpha.rawptr(), tss->T, J21.rawptr());
  J21 *= dg;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,n)] = -J21.rawptr()[CINDEX(i,j,6)];
    }
  }

  // J22
  History J22 = gather_blank_state_().derivative<History>();
  flow_->dh_da(s_np1.data(), alpha.rawptr(), tss->T, J22.rawptr());
  J22 *= dg;
  for (int i=0; i<nh; i++) J22.rawptr()[CINDEX(i,i,nh)] -= 1.0;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<nh; j++) {
      J[CINDEX((i+6), (j+6), n)] = -J22.rawptr()[CINDEX(i,j,nh)];
    }
  }

  // J23
  for (int i=0; i<nh; i++) {
    J[CINDEX((i+6), (6+nh), n)] = -h.rawptr()[i];
  }

  // J31
  Symmetric J31;
  flow_->df_ds(s_np1.data(), alpha.rawptr(), tss->T, J31.s());
  for (int i=0; i<6; i++) {
    J[CINDEX((6+nh), i, n)] = J31(i);
  }

  // J32
  History J32 = gather_blank_state_();
  flow_->df_da(s_np1.data(), alpha.rawptr(), tss->T, J32.rawptr());
  for (int i=0; i<nh; i++) {
    J[CINDEX((6+nh), (i+6), n)] = J32.rawptr()[i];
  }

  // J33
  J[CINDEX((6+nh), (6+nh), n)] = 0.0;
}

const std::shared_ptr<const LinearElasticModel> SmallStrainRateIndependentPlasticity::elastic() const
{
  return elastic_;
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
  hist.add<Symmetric>(prefix("plastic_strain"));
  plastic_->set_variable_prefix(get_variable_prefix());
  plastic_->populate_hist(hist);
}

void SmallStrainCreepPlasticity::init_state(History & hist) const
{
  hist.get<Symmetric>(prefix("plastic_strain")) = Symmetric::zero();
  plastic_->init_hist(hist);
}

void SmallStrainCreepPlasticity::update_sd_state(
    const Symmetric & E_np1, const Symmetric & E_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    Symmetric & S_np1, const Symmetric & S_n,
    History & H_np1, const History & H_n,
    SymSymR4 & AA_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Split out the history
  History base_np1 = plastic_->gather_state_(H_np1.rawptr());
  History base_n = plastic_->gather_state_(H_n.rawptr());
  
  // Solve the system to get the update
  std::unique_ptr<SSCPTrialState> ts = make_trial_state(E_np1, E_n, 
                                       T_np1, T_n, t_np1, t_n, S_n,
                                       H_n);

  std::vector<double> xv(nparams());
  double * x = &xv[0];
  solve(this, x, ts.get(), {rtol_, atol_, miter_, verbose_, linesearch_});

  // Store the ep strain
  H_np1.get<Symmetric>(prefix("plastic_strain")) = Symmetric(x);

  // Do the plastic update to get the new history and stress
  SymSymR4 A;
  plastic_->update_sd_state(Symmetric(x), ts->ep_strain,
                            T_np1, T_n,
                            t_np1, t_n, S_np1, S_n,
                            base_np1, base_n,
                            A, u_np1, u_n, p_np1, p_n);

  // Do the creep update to get a tangent component
  Symmetric creep_old = E_n - ts->ep_strain;
  Symmetric creep_new;
  SymSymR4 B;
  creep_->update(S_np1.data(), creep_new.s(), creep_old.s(), T_np1, T_n,
                 t_np1, t_n, B.s());

  // Form the relatively simple tangent
  AA_np1 = form_tangent_(A, B);

  // Energy calculation (trapezoid rule)
  auto [du, dp] = trapezoid_energy(E_np1, E_n, creep_new, creep_old, S_np1, S_n);
  u_np1 = u_n + du;
  p_np1 += dp;
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
  std::copy(tss->ep_strain.data(), tss->ep_strain.data() + 6, x);
}

void SmallStrainCreepPlasticity::RJ(const double * const x, TrialState * ts, 
                                   double * const R, double * const J)
{
  SSCPTrialState * tss = static_cast<SSCPTrialState*>(ts);

  // First update the elastic-plastic model
  Symmetric S_np1;
  SymSymR4 A_np1;
  History H_np1 = plastic_->gather_blank_state_();

  double u_np1, u_n;
  double p_np1, p_n;
  u_n = 0.0;
  p_n = 0.0;

  Symmetric ep_np1(x);
  plastic_->update_sd_state(ep_np1, tss->ep_strain,
                            tss->T_np1, tss->T_n,
                            tss->t_np1, tss->t_n, S_np1, 
                            tss->s_n,
                            H_np1, tss->h_n, A_np1,
                            u_np1, u_n, p_np1, p_n);

  // Then update the creep strain
  Symmetric creep_old = tss->e_n - tss->ep_strain;
  Symmetric creep_new;
  SymSymR4 B;
  creep_->update(S_np1.data(), creep_new.s(), creep_old.data(), tss->T_np1, tss->T_n,
                 tss->t_np1, tss->t_n, B.s());

  // Form the residual
  Symmetric Rs(R);
  Rs = (ep_np1 + creep_new - tss->e_np1) * sf_;
  
  // The Jacobian is a straightforward combination of the two derivatives
  SymSymR4 Js(J);
  Js = (B.dot(A_np1) + SymSymR4::id()) * sf_;
}

std::unique_ptr<SSCPTrialState>
SmallStrainCreepPlasticity::make_trial_state(
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const Symmetric & s_n, const History & h_n)
{
  return std::make_unique<SSCPTrialState>(
      h_n.get<Symmetric>(prefix("plastic_strain")),
      e_n, e_np1, s_n, T_n, T_np1, t_n, t_np1,
      plastic_->gather_state_(h_n.rawptr()));
}

SymSymR4 SmallStrainCreepPlasticity::form_tangent_(
    const SymSymR4 & A, const SymSymR4 & B)
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
  return A - A.dot(B.dot((A.dot(B) + SymSymR4::id()).inverse().dot(A)));
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

std::unique_ptr<TrialState> GeneralIntegrator::setup(
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const Symmetric & s_n,
    const History & h_n)
{
  // Time increment
  double dt = t_np1 - t_n;
  
  // Rates, avoiding divide by zero
  double T_dot = 0.0;
  Symmetric e_dot;
  if (dt > 0.0) {
    T_dot = (T_np1 - T_n) / dt;
    e_dot = (e_np1 - e_n) / dt;
  }
  
  // Default to elastic guess, unless the skip_first_ flag is true
  Symmetric s_guess;
  if ((t_n > 0.0) || (!skip_first_)) 
    s_guess = elastic_->C(T_np1).dot(e_np1 - e_n) + s_n;
  
  return std::make_unique<GITrialState>(e_dot, s_n, s_guess, h_n, T_np1, 
                                        T_dot, dt);
}

bool GeneralIntegrator::elastic_step(
    const TrialState * ts, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n, const Symmetric & s_n,
    const History & h_n)
{
  double dt = t_np1 - t_n;

  if (dt < std::numeric_limits<double>::epsilon()) {
    return true;
  }

  return false;
}

void GeneralIntegrator::update_internal(
    const double * const x, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    Symmetric & s_np1, const Symmetric & s_n,
    History & h_np1, const History & h_n)
{
  s_np1.copy_data(x);
  h_np1.copy_data(x+6);
}

void GeneralIntegrator::strain_partial(
    const TrialState * ts, const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const Symmetric & s_np1, const Symmetric & s_n, 
    const History & h_np1, const History & h_n, double * de)
{
  const GITrialState * tss = static_cast<const GITrialState*>(ts);
  
  SymSymR4 estress;
  rule_->ds_de(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(), tss->T, tss->Tdot, estress.s());
  
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      de[CINDEX(i,j,6)] = estress(i,j);
    }
  }
  
  History ehist = gather_blank_state_().derivative<Symmetric>();
  rule_->da_de(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(), tss->T, tss->Tdot, ehist.rawptr());
  for (size_t i = 0; i < nstate(); i++) {
    for (size_t j = 0; j < 6; j++) {
      de[CINDEX((i+6),j,6)] = ehist.rawptr()[CINDEX(i,j,6)];
    }
  }
}

void GeneralIntegrator::work_and_energy(
    const TrialState * ts,
    const Symmetric & e_np1, const Symmetric & e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    const Symmetric & s_np1, const Symmetric & s_n,
    const History & h_np1, const History & h_n,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  const GITrialState * tss = static_cast<const GITrialState *>(ts);

  // Energy
  Symmetric de = e_np1 - e_n;
  Symmetric ds = s_np1 - s_n;
  u_np1 = u_n + ds.contract(de) / 2.0;

  // Dissipation needs a special call
  double p_dot_np1;
  rule_->work_rate(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(),
                   T_np1, tss->Tdot, p_dot_np1);
  double p_dot_n;
  rule_->work_rate(s_n.data(), h_n.rawptr(), 
                   tss->e_dot.data(), T_n, tss->Tdot, p_dot_n);
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
  std::copy(tss->s_guess.data(), tss->s_guess.data()+6, x);
  std::copy(tss->h_n.rawptr(), tss->h_n.rawptr() + nstate(), &x[6]);

  rule_->override_guess(x);
}

void GeneralIntegrator::RJ(const double * const x, TrialState * ts,
                          double * const R, double * const J)
{
  GITrialState * tss = static_cast<GITrialState*>(ts);

  // Setup
  Symmetric s_np1(x);
  History h_np1 = gather_state_(x+6);
  
  // Helps with vectorization
  // Really as I declared both const this shouldn't be necessary but hey
  // I don't design optimizing compilers for a living
  int nstate = this->nstate();
  int nparams = this->nparams();

  // Residual calculation
  Symmetric fr;
  rule_->s(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(), tss->T, tss->Tdot, fr.s());
  Symmetric R1 = s_np1 - tss->s_n - fr * tss->dt;
  for (int i=0; i<6; i++) {
    R[i] = R1(i);
  }
  History h_rate = gather_blank_state_();
  rule_->a(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(), tss->T, tss->Tdot, h_rate.rawptr());
  for (int i=0; i<nstate; i++) {
    R[i+6] = h_np1.rawptr()[i] - tss->h_n.rawptr()[i] - h_rate.rawptr()[i] * tss->dt;
  }

  // Jacobian calculation
  SymSymR4 J11;
  rule_->ds_ds(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(), tss->T, tss->Tdot, J11.s());
  J11 = J11 * tss->dt - SymSymR4::id();
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,nparams)] = -J11.data()[CINDEX(i,j,6)];
    }
  }
  
  // Transpose?
  History J12 = gather_blank_state_().derivative<Symmetric>();
  rule_->ds_da(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(), tss->T, tss->Tdot, J12.rawptr());
  for (int i=0; i<6; i++) {
    for (int j=0; j<nstate; j++) {
      J[CINDEX(i,(j+6),nparams)] = -J12.rawptr()[CINDEX(i,j,nstate)] * tss->dt;
    }
  }
  
  // Transpose?
  History J21 = gather_blank_state_().derivative<Symmetric>();
  rule_->da_ds(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(), tss->T, tss->Tdot, J21.rawptr());
  for (int i=0; i<nstate; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,nparams)] = -J21.rawptr()[CINDEX(i,j,6)] * tss->dt;
    }
  }
  
  History J22 = gather_blank_state_().derivative<History>();
  rule_->da_da(s_np1.data(), h_np1.rawptr(), tss->e_dot.data(), tss->T, tss->Tdot, J22.rawptr());

  // More vectorization
  double dt = tss->dt;
  for (int i=0; i<nstate*nstate; i++) J22.rawptr()[i] *= dt;
  for (int i=0; i<nstate; i++) J22.rawptr()[CINDEX(i,i,nstate)] -= 1.0;

  for (int i=0; i<nstate; i++) {
    for (int j=0; j<nstate; j++) {
      J[CINDEX((i+6),(j+6),nparams)] = -J22.rawptr()[CINDEX(i,j,nstate)];
    }
  }
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

void KMRegimeModel::update_sd_state(
    const Symmetric & E_np1, const Symmetric & E_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    Symmetric & S_np1, const Symmetric & S_n,
    History & H_np1, const History & H_n,
    SymSymR4 & AA_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  // Calculate activation energy
  double g = activation_energy_(E_np1, E_n, T_np1, t_np1, t_n);

  // Note this relies on everything being sorted.  You probably want to
  // error check at some point
  for (size_t i=0; i<gs_.size(); i++) {
    if (g < gs_[i]) {
      return models_[i]->update_sd_state(E_np1, E_n, T_np1, T_n, t_np1, t_n, 
                                   S_np1, S_n, H_np1, H_n, AA_np1, u_np1, u_n,
                                   p_np1, p_n);
    }
  }
  return models_.back()->update_sd_state(E_np1, E_n, T_np1, T_n, t_np1, t_n,
                                   S_np1, S_n, H_np1, H_n, AA_np1, u_np1, u_n,
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

double KMRegimeModel::activation_energy_(const Symmetric & e_np1, 
                                         const Symmetric & e_n,
                                         double T_np1,
                                         double t_np1, double t_n)
{
  Symmetric de = (e_np1 - e_n) / (t_np1 - t_n);
  double rate = sqrt(2.0/3.0) * de.norm();
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

std::tuple<double,double> trapezoid_energy(
    const Symmetric & e_np1, const Symmetric & e_n,
    const Symmetric & ep_np1, const Symmetric & ep_n,
    const Symmetric & s_np1, const Symmetric & s_n)
{
  Symmetric de = e_np1 - e_n;
  Symmetric dp = ep_np1 - ep_n;
  Symmetric ds = s_np1 + s_n;
  
  double du = ds.contract(de)/2.0;
  double dw = ds.contract(dp)/2.0;

  return std::tie(du,dw);
}

} // namespace neml
