#include "neml.h"

#include "nemlmath.h"
#include "nemlerror.h"

#include <cassert>

#include <iostream>

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
    double * const A_np1)
{
  assert(false); // Pass on implementing for now
}

int NEMLModel_ldF::update_sd(
     const double * const e_np1, const double * const e_n,
     double T_np1, double T_n,
     double t_np1, double t_n,
     double * const s_np1, const double * const s_n,
     double * const h_np1, const double * const h_n,
     double * const A_np1)
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
    double * const A_np1)
{
  assert(false); // Pass for now
}

int NEMLModel_ldI::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1)
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
    double * const A_np1)
{
  assert(false); // Pass for now
}

int NEMLModel_sd::update_ldI(
    const double * const l_inc,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1)
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
       double * const A_np1)
{
  elastic_->C(T_np1, A_np1);
  mat_vec(A_np1, 6, e_np1, 6, s_np1);
  return 0;
}


// Implementation of small strain rate independent plasticity
//

SmallStrainRateIndependentPlasticity::SmallStrainRateIndependentPlasticity(
    std::shared_ptr<LinearElasticModel> elastic,
    std::shared_ptr<RateIndependentFlowRule> flow, double rtol, double atol,
    int miter, bool verbose) :
      elastic_(elastic), flow_(flow), rtol_(rtol), atol_(atol), miter_(miter),
      verbose_(verbose)
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
       double * const A_np1)
{
  // Setup and store the trial state for the solver
  set_trial_state(e_np1, h_n, T_np1);

  // Check to see if this is an elastic state
  double fv;
  flow_->f(s_tr_, &h_tr_[0], T_np1, fv);
  double dg;

  // If elastic, copy over and return
  if (fv < atol_) {
    std::copy(s_tr_, s_tr_+6, s_np1);
    std::copy(ep_tr_, ep_tr_+6, h_np1);
    std::copy(&h_tr_[0], &h_tr_[0]+flow_->nhist(), &h_np1[6]);
    std::copy(C_, C_+36, A_np1);
    dg = 0.0;
    return 0;
  }
  // Else solve and extract updated parameters from the solver vector
  else {
    double x[nparams()];
    int ier = solve(shared_from_this(), x, atol_, rtol_, miter_, verbose_);
    if (ier != SUCCESS) return ier;

    // Extract solved parameters
    std::copy(x, x+6, h_np1); // plastic strain
    std::copy(x+6, x+6+flow_->nhist(), &h_np1[6]); // history
    dg = x[6+flow_->nhist()];
    double ee[6];
    sub_vec(e_np1, h_np1, 6, ee);
    mat_vec(C_, 6, ee, 6, s_np1);

    // Complicated tangent calc...
    
    // Check K-T and return
    return 0;
  }
}

size_t SmallStrainRateIndependentPlasticity::nparams() const
{
  // My convention: plastic strain, history, consistency parameter
  return 6 + flow_->nhist() + 1;
}

int SmallStrainRateIndependentPlasticity::init_x(double * const x)
{
  std::copy(ep_tr_, ep_tr_+6, x);
  std::copy(&h_tr_[0], &h_tr_[0]+flow_->nhist(), &x[6]);
  x[6+flow_->nhist()] = 0.0; // consistency parameter
  return 0;
}

int SmallStrainRateIndependentPlasticity::RJ(const double * const x, double * const R, double * const J)
{
  // Setup from current state
  const double * const ep = &x[0];
  const double * const alpha  = &x[6];
  const double & dg = x[6+flow_->nhist()];
  double ee[6];
  sub_vec(e_np1_, ep, 6, ee);
  double s[6];
  mat_vec(C_, 6, ee, 6, s);

  // Residual calculation
  double g[6];
  flow_->g(s, alpha, T_, g); 
  double h[flow_->nhist()];
  flow_->h(s, alpha, T_, h);
  double f;
  flow_->f(s, alpha, T_, f);
  for (int i=0; i<6; i++) 
  {
    R[i] = -ep[i] + ep_tr_[i] + g[i] * dg;
  }
  for (int i=0; i<flow_->nhist(); i++) {
    R[i+6] = -alpha[i] + h_tr_[i] + h[i] * dg;
  }
  R[6+flow_->nhist()] = f;

  // Now the jacobian calculation...
  int n = nparams();
  int nh = flow_->nhist();
  
  // J11
  double gs[36];
  double J11[36];
  flow_->dg_ds(s, alpha, T_, gs);
  mat_mat(6, 6, 6, gs, C_, J11);
  for (int i=0; i<36; i++) J11[i] = -J11[i] * dg;
  for (int i=0; i<6; i++) J11[CINDEX(i,i,6)] -= 1.0;
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,n)] = J11[CINDEX(i,j,6)];
    }
  }
  
  // J12
  double J12[6*nh];
  flow_->dg_da(s, alpha, T_, J12);
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
  flow_->dh_ds(s, alpha, T_, ha);
  mat_mat(nh, 6, 6, ha, C_, J21);
  for (int i=0; i<nh*6; i++) J21[i] = -J21[i] * dg;
  for (int i=0; i<nh; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX((i+6),j,n)] = J21[CINDEX(i,j,6)];
    }
  }

  // J22
  double J22[nh*nh];
  flow_->dh_da(s, alpha, T_, J22);
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
  flow_->df_ds(s, alpha, T_, fs);
  mat_vec_trans(C_, 6, fs, 6, J31);
  for (int i=0; i<6; i++) J31[i] = -J31[i];
  for (int i=0; i<6; i++) {
    J[CINDEX((6+nh), i, n)] = J31[i];
  }

  // J32
  double J32[nh];
  flow_->df_da(s, alpha, T_, J32);
  for (int i=0; i<nh; i++) {
    J[CINDEX((6+nh), (i+6), n)] = J32[i];
  }

  // J33
  J[CINDEX((6+nh), (6+nh), n)] = 0.0;
  
  return 0;
}

int SmallStrainRateIndependentPlasticity::set_trial_state(
    const double * const e_np1, const double * const h_n, double T_np1)
{
  // Save e_np1
  std::copy(e_np1, e_np1+6, e_np1_);
  // ep_tr = ep_n
  std::copy(h_n, h_n+6, ep_tr_);
  // h_tr = h_n
  h_tr_.resize(flow_->nhist());
  std::copy(h_n+6, h_n+nhist(), h_tr_.begin());
  // Calculate the trial stress
  double ee[6];
  sub_vec(e_np1, ep_tr_, 6, ee);
  elastic_->C(T_np1, C_);
  mat_vec(C_, 6, ee, 6, s_tr_);
  // Store temp
  T_ = T_np1;
  return 0;
}


} // namespace neml
