#include "neml.h"

#include "nemlmath.h"
#include "nemlerror.h"

#include <cassert>
#include <limits>

namespace neml {

// Base class NEMLModel

NEMLModel::NEMLModel()
{

}

NEMLModel::~NEMLModel()
{

}



// NEMLModel_ldF implementation
NEMLModel_ldF::NEMLModel_ldF() :
    NEMLModel()
{

}

NEMLModel_ldF::~NEMLModel_ldF()
{

}

size_t NEMLModel_ldF::nstore() const
{
  return nhist() + 0;
}

int NEMLModel_ldF::init_store(double * const store) const
{
  assert(false); 
}

int NEMLModel_ldF::update_ldI(
    const double * const l_inc,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  assert(false); // Pass on implementing for now
}

int NEMLModel_ldF::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  assert(false); // Pass on implementing for now
}

// NEMLModel_ldI implementation
NEMLModel_ldI::NEMLModel_ldI() :
    NEMLModel()
{

}

NEMLModel_ldI::~NEMLModel_ldI()
{

}

size_t NEMLModel_ldI::nstore() const
{
  return nhist() + 0;
}

int NEMLModel_ldI::init_store(double * const store) const
{
  assert(false); 
}

int NEMLModel_ldI::update_ldF(
    const double * const F_np1, const double * const F_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  assert(false); // Pass for now
}

int NEMLModel_ldI::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  assert(false); // Pass for now
}


// NEMLModel_sd implementation
NEMLModel_sd::NEMLModel_sd() :
    NEMLModel()
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
}

int NEMLModel_sd::update_ldF(
    const double * const F_np1, const double * const F_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  assert(false); // Pass for now
}

int NEMLModel_sd::update_ldI(
    const double * const l_inc,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  assert(false); // Pass for now
}


// Implementation of small strain elasticity
SmallStrainElasticity::SmallStrainElasticity(std::shared_ptr<LinearElasticModel> elastic) :
    elastic_(elastic)
{


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

  // Energy calculation
  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  u_np1 = u_n + dot_vec(s_np1, de, 6);
  p_np1 = p_n;

  return 0;
}


// Implementation of small strain rate independent plasticity
//

SmallStrainRateIndependentPlasticity::SmallStrainRateIndependentPlasticity(
    std::shared_ptr<LinearElasticModel> elastic,
    std::shared_ptr<RateIndependentFlowRule> flow, double tol,
    int miter, bool verbose, double kttol, bool check_kt) :
      elastic_(elastic), flow_(flow), tol_(tol), miter_(miter),
      verbose_(verbose), kttol_(kttol), check_kt_(check_kt)
{

}

size_t SmallStrainRateIndependentPlasticity::nhist() const
{
  // My convention: plastic strain then actual history
  return 6 + flow_->nhist();
}

int SmallStrainRateIndependentPlasticity::init_hist(double * const hist) const
{
  std::fill(hist, hist+6, 0.0);
  return flow_->init_hist(&hist[6]);
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
  make_trial_state(e_np1, h_n, T_np1, t_np1, t_n, ts);

  // Check to see if this is an elastic state
  double fv;
  flow_->f(ts.s_tr, &ts.h_tr[0], T_np1, fv);
  double dg;

  // If elastic, copy over and return
  if (fv < tol_) {
    std::copy(ts.s_tr, ts.s_tr+6, s_np1);
    std::copy(ts.ep_tr, ts.ep_tr+6, h_np1);
    std::copy(&ts.h_tr[0], &ts.h_tr[0]+flow_->nhist(), &h_np1[6]);
    std::copy(ts.C, ts.C+36, A_np1);
    dg = 0.0;
  }
  // Else solve and extract updated parameters from the solver vector
  else {
    double x[nparams()];
    int ier = solve(this, x, &ts, tol_, miter_, verbose_);
    if (ier != SUCCESS) return ier;

    // Extract solved parameters
    std::copy(x, x+6, h_np1); // plastic strain
    std::copy(x+6, x+6+flow_->nhist(), &h_np1[6]); // history
    dg = x[6+flow_->nhist()];
    double ee[6];
    sub_vec(e_np1, h_np1, 6, ee);
    mat_vec(ts.C, 6, ee, 6, s_np1);

    // Complicated tangent calc...
    calc_tangent_(x, &ts, s_np1, h_np1, dg, A_np1);
  }

  // Energy calculation
  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  u_np1 = u_n + dot_vec(s_np1, de, 6);

  double dep[6];
  sub_vec(h_np1, h_n, 6, dep);
  p_np1 = p_n + dot_vec(s_np1, dep, 6);

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
  double h[flow_->nhist()];
  flow_->h(s, alpha, tss->T, h);
  double f;
  flow_->f(s, alpha, tss->T, f);
  for (int i=0; i<6; i++) 
  {
    R[i] = -ep[i] + tss->ep_tr[i] + g[i] * dg;
  }
  for (int i=0; i<flow_->nhist(); i++) {
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
  double J12[6*nh];
  flow_->dg_da(s, alpha, tss->T, J12);
  for (int i=0; i<6*nh; i++) J12[i] *= dg;
  for (int i=0; i<6; i++) {
    for (int j=0; j < nh; j++) {
      J[CINDEX(i,(j+6),n)] = J12[CINDEX(i,j,nh)];
    }
  }

  // J13
  for (int i=0; i<6; i++) {
    J[CINDEX(i,(6+nh),n)] = g[i];
  }

  // J21
  double ha[nh*6];
  double J21[nh*6];
  flow_->dh_ds(s, alpha, tss->T, ha);
  mat_mat(nh, 6, 6, ha, tss->C, J21);
  for (int i=0; i<nh*6; i++) J21[i] = -J21[i] * dg;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,n)] = J21[CINDEX(i,j,6)];
    }
  }

  // J22
  double J22[nh*nh];
  flow_->dh_da(s, alpha, tss->T, J22);
  for (int i=0; i<nh*nh; i++) J22[i] *= dg;
  for (int i=0; i<nh; i++) J22[CINDEX(i,i,nh)] -= 1.0;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<nh; j++) {
      J[CINDEX((i+6), (j+6), n)] = J22[CINDEX(i,j,nh)];
    }
  }

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
  double J32[nh];
  flow_->df_da(s, alpha, tss->T, J32);
  for (int i=0; i<nh; i++) {
    J[CINDEX((6+nh), (i+6), n)] = J32[i];
  }

  // J33
  J[CINDEX((6+nh), (6+nh), n)] = 0.0;
  
  return 0;
}

int SmallStrainRateIndependentPlasticity::make_trial_state(
    const double * const e_np1, const double * const h_n, double T_np1,
    double t_np1, double t_n, SSRIPTrialState & ts)
{
  // Save e_np1
  std::copy(e_np1, e_np1+6, ts.e_np1);
  // ep_tr = ep_n
  std::copy(h_n, h_n+6, ts.ep_tr);
  // h_tr = h_n
  ts.h_tr.resize(flow_->nhist());
  std::copy(h_n+6, h_n+nhist(), ts.h_tr.begin());
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

  double R[nparams()];
  double J[nparams()*nparams()];
  
  RJ(x, ts, R, J);

  int n = nparams();
  int nk = 6;
  int ne = nparams() - nk;
  
  double Jkk[nk*nk];
  for (int i=0; i<nk; i++) {
    for (int j=0; j<nk; j++) {
      Jkk[CINDEX(i,j,nk)] = J[CINDEX(i,j,n)];
    }
  }

  double Jke[nk*ne];
  for (int i=0; i<nk; i++) {
    for (int j=0; j<ne; j++) {
      Jke[CINDEX(i,j,ne)] = J[CINDEX(i,(j+nk),n)];
    }
  }

  double Jek[ne*nk];
  for (int i=0; i<ne; i++) {
    for (int j=0; j<nk; j++) {
      Jek[CINDEX(i,j,nk)] = J[CINDEX((i+nk),(j),n)];
    }
  }

  double Jee[ne*ne];
  for (int i=0; i<ne; i++) {
    for (int j=0; j<ne; j++) {
      Jee[CINDEX(i,j,ne)] = J[CINDEX((i+nk),(j+nk),n)];
    }
  }

  invert_mat(Jee, ne);
  
  double A[nk*6];
  double B[ne*6];
  
  int nh = flow_->nhist();

  double dg_ds[6*6];
  flow_->dg_ds(s_np1, &h_np1[6], tss->T, dg_ds);
  double dh_ds[nh*6];
  flow_->dh_ds(s_np1, &h_np1[6], tss->T, dh_ds);
  double df_ds[6];
  flow_->df_ds(s_np1, &h_np1[6], tss->T, df_ds);

  mat_mat(6, 6, 6, dg_ds, tss->C, A);
  for (int i=0; i<nk*6; i++) A[i] *= dg;
  
  mat_mat(nh, 6, 6, dh_ds, tss->C, B);
  for (int i=0; i<nh*6; i++) B[i] *= dg;
  mat_vec_trans(tss->C, 6, df_ds, 6, &B[nh*6]);
  
  double T1[ne*nk];
  mat_mat(ne, nk, ne, Jee, Jek, T1);
  double T2[nk*nk];
  mat_mat(nk, nk, ne, Jke, T1, T2);
  for (int i=0; i<nk*nk; i++) T2[i] = Jkk[i] - T2[i];
  invert_mat(T2, nk);

  double T3[ne*nk];
  mat_mat(ne, nk, ne, Jee, B, T3);
  double T4[nk*nk];
  mat_mat(nk, nk, ne, Jke, T3, T4);
  for (int i=0; i<nk*nk; i++) T4[i] -= A[i];

  double dep[36];
  mat_mat(6, 6, 6, T2, T4, dep);

  for (int i=0; i<36; i++) dep[i] = -dep[i];
  for (int i=0; i<6; i++) dep[CINDEX(i,i,6)] += 1.0;

  mat_mat(6, 6, 6, tss->C, dep, A_np1);

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

// Start general integrator implementation
GeneralIntegrator::GeneralIntegrator(std::shared_ptr<GeneralFlowRule> rule,
                                     double tol, int miter,
                                     bool verbose, int max_divide) :
    rule_(rule), tol_(tol), miter_(miter), verbose_(verbose), 
    max_divide_(max_divide)
{

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
    make_trial_state(e_next, e_past, s_past, h_past, T_next, T_past, t_next, 
                     t_past, ts);

    // Solve for x
    double x[nparams()];
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
  make_trial_state(e_np1, e_n, s_n, h_n, T_np1, T_n, t_np1, t_n, ts);
  double y[nparams()];
  std::copy(s_np1, s_np1+6, y);
  std::copy(h_np1, h_np1+nhist(), &y[6]);
  
  calc_tangent_(y, &ts, A_np1);
  
  // Energy calculation
  double de[6];
  sub_vec(e_np1, e_n, 6, de);
  u_np1 = u_n + dot_vec(s_np1, de, 6);

  // FIXME
  p_np1 = p_n;

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
  double s_mod[36];
  std::copy(x, x+6, s_mod);
  if (norm2_vec(x, 6) < std::numeric_limits<double>::epsilon()) {
    s_mod[0] = 2.0 * std::numeric_limits<double>::epsilon();
  }
  const double * const h_np1 = &x[6];


  // Residual calculation
  rule_->s(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, R);
  for (int i=0; i<6; i++) {
    R[i] = -s_mod[i] + tss->s_n[i] + R[i] * tss->dt;
  }
  rule_->a(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, &R[6]);
  for (int i=0; i<nhist(); i++) {
    R[i+6] = -h_np1[i] + tss->h_n[i] + R[i+6] * tss->dt;
  }

  // Jacobian calculation
  double J11[36];
  rule_->ds_ds(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, J11);
  for (int i=0; i<36; i++) J11[i] *= tss->dt;
  for (int i=0; i<6; i++) J11[CINDEX(i,i,6)] -= 1.0;
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,nparams())] = J11[CINDEX(i,j,6)];
    }
  }

  double J12[6*nhist()];
  rule_->ds_da(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, J12);
  for (int i=0; i<6; i++) {
    for (int j=0; j<nhist(); j++) {
      J[CINDEX(i,(j+6),nparams())] = J12[CINDEX(i,j,nhist())] * tss->dt;
    }
  }

  double J21[nhist()*6];
  rule_->da_ds(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, J21);
  for (int i=0; i<nhist(); i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,nparams())] = J21[CINDEX(i,j,6)] * tss->dt;
    }
  }

  double J22[nhist()*nhist()];
  rule_->da_da(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, J22);
  for (int i=0; i<nhist()*nhist(); i++) J22[i] *= tss->dt;
  for (int i=0; i<nhist(); i++) J22[CINDEX(i,i,nhist())] -= 1.0;

  for (int i=0; i<nhist(); i++) {
    for (int j=0; j<nhist(); j++) {
      J[CINDEX((i+6),(j+6),nparams())] = J22[CINDEX(i,j,nhist())];
    }
  }

  return 0;
}


int GeneralIntegrator::make_trial_state(
    const double * const e_np1, const double * const e_n,
    const double * const s_n, const double * const h_n,
    double T_np1, double T_n, double t_np1, double t_n,
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

  // Call for extra derivatives
  double A[36];
  rule_->ds_de(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, A);
  double B[nhist()*6];
  rule_->da_de(s_mod, h_np1, tss->e_dot, tss->T, tss->Tdot, B);

  // Call for the jacobian
  double R[nparams()];
  double J[nparams()*nparams()];
  RJ(x, ts, R, J);

  // Separate blocks...
  int n = nparams();
  
  double J11[6*6];
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J11[CINDEX(i,j,6)] = J[CINDEX(i,j,n)];
    }
  }

  double J12[6*nhist()];
  for (int i=0; i<6; i++) {
    for (int j=0; j<nhist(); j++) {
      J12[CINDEX(i,j,nhist())] = J[CINDEX(i,(j+6),n)];
    }
  }

  double J21[nhist()*6];
  for (int i=0; i<nhist(); i++) {
    for (int j=0; j<6; j++) {
      J21[CINDEX(i,j,6)] = J[CINDEX((i+6),(j),n)];
    }
  }

  double J22[nhist()*nhist()];
  for (int i=0; i<nhist(); i++) {
    for (int j=0; j<nhist(); j++) {
      J22[CINDEX(i,j,nhist())] = J[CINDEX((i+6),(j+6),n)];
    }
  }

  // Only need the inverse
  invert_mat(J22, nhist());

  // Start multiplying through
  double T1[nhist()*6];
  mat_mat(nhist(), 6, nhist(), J22, J21, T1);
  double T2[6*6];
  mat_mat(6, 6, nhist(), J12, T1, T2);
  for (int i=0; i<6*6; i++) T2[i] = J11[i] - T2[i];
  invert_mat(T2, 6);

  double T3[nhist()*6];
  mat_mat(nhist(), 6, nhist(), J22, B, T3);
  double T4[6*6];
  mat_mat(6, 6, nhist(), J12, T3, T4);
  for (int i=0; i<6*6; i++) T4[i] -= A[i];

  mat_mat(6, 6, 6, T2, T4, A_np1);

  return 0;
}


} // namespace nhist()ml
