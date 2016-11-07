#include "deviatoric.h"

#include "nemlmath.h"
#include "nemlerror.h"

#include <limits>
#include <algorithm>
#include <iostream>

namespace neml {

// Base class implementation
DeviatoricModel::DeviatoricModel()
{

}

DeviatoricModel::~DeviatoricModel()
{

}

int DeviatoricModel::get_dev(const double * const in, double * const out) const
{
  double mean = in[0] + in[1] + in[2];
  std::copy(in, in+6, out);
  for (int i=0; i<3; i++) {
    out[i] -= mean / 3.0;
  }

  return 0;
}

// Linear elastic model
LEModel::LEModel(std::shared_ptr<ShearModulus> modulus) :
    modulus_(modulus)
{

}

LEModel::~LEModel()
{

}

size_t LEModel::nhist() const
{
  return 0;
}

int LEModel::init_hist(double* const h) const
{
  return 0;
}

int LEModel::update(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{
  double edev[6];
  get_dev(e_np1, edev);
  std::copy(edev, edev+6, s_np1);
  double mu = modulus_->modulus(T_np1);
  for (int i=0; i<6; i++) {
    s_np1[i] *= 2.0 * mu;
  }

  std::fill(A_np1,A_np1+36,0.0);

  A_np1[CINDEX(0,0,6)] = 4.0 / 3.0 * mu;
  A_np1[CINDEX(0,1,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(0,2,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(1,0,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(1,1,6)] = 4.0 / 3.0 * mu;
  A_np1[CINDEX(1,2,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(2,0,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(2,1,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(2,2,6)] = 4.0 / 3.0 * mu;
  A_np1[CINDEX(3,3,6)] = 2.0 * mu;
  A_np1[CINDEX(4,4,6)] = 2.0 * mu;
  A_np1[CINDEX(5,5,6)] = 2.0 * mu;

  return 0;
}

// Rate independent, associative flow model implementation
RIAFModel::RIAFModel(std::shared_ptr<ShearModulus> modulus,
                     std::shared_ptr<YieldSurface> surface, 
                     std::shared_ptr<AssociativeHardening> hardening,
                     double tol, int miter, bool verbose) :
    modulus_(modulus), surface_(surface), hardening_(hardening),
    tol_(tol), miter_(miter), verbose_(verbose)
{
}

RIAFModel::~RIAFModel()
{

}

size_t RIAFModel::nhist() const
{

  return hardening_->nhist() + 6; // Need to keep around the plastic strain
}

int RIAFModel::init_hist(double* const h) const
{
  if (hardening_->nhist() != surface_->nhist()) {
    return INCOMPATIBLE_MODELS;
  }

  int ier = hardening_->init_hist(h);
  if (ier != 0) {
    return ier;
  }

  std::fill(h+hardening_->nhist(), h+hardening_->nhist()+6, 0.0);

  return 0;
}


int RIAFModel::update(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{
  // Get setup in deviatoric space
  double e_dev_np1[6];
  get_dev(e_np1, e_dev_np1);

  // Setup iterators
  std::copy(h_n, h_n+nhist(), h_np1);
  double * const alpha = &h_np1[0];
  double * const ep = &h_np1[hardening_->nhist()];

  const double * const alpha_n = &h_n[0];
  const double * const ep_n = &h_n[hardening_->nhist()];

  double dg = 0;

  // Shear modulus will remain constant over the step
  double mu = modulus_->modulus(T_np1);

  // Iterate until convergence
  int i;
  if (verbose_) {
    std::cout << "Iter.\tNorm.\tYield" << std::endl;
  }
  for (i=0; i<miter_; i++) {
    if (verbose_) std::cout << i;
    if (take_step(e_dev_np1, mu, T_np1, ep, alpha, ep_n, alpha_n, dg, 
                  s_np1)) break;
  }
  
  if (i == miter_) return MAX_ITERATIONS;
  
  // Setup the algorithmic tangent
  get_tangent(mu, T_np1, alpha, dg, s_np1, A_np1);

  return 0;
}

bool RIAFModel::take_step(const double * const e_dev, double mu, double T,
               double * const ep_np1, double * const alpha_np1,
               const double * const ep_n, const double * const alpha_n,
               double & dg, double * const s) const
{
  // Setup
  for (int i=0; i<6; i++) {
    s[i] = 2.0 * mu * (e_dev[i] - ep_np1[i]);
  }
  
  int nh = hardening_->nhist();
  double q[nh];
  hardening_->q(alpha_np1, T, q);

  double f;
  surface_->f(s, q, T, f);
  
  // Yield surface gradient vector
  double G[nhist()];
  surface_->df_ds(s, q, T, G);
  surface_->df_dq(s, q, T, &G[6]);

  // Inverse of plastic moduli
  double Di[nh*nh];
  hardening_->D_inv(alpha_np1, T, Di);

  // Residual
  double R[nhist()];
  for (int i=0; i<6; i++) {
    R[i] = -ep_np1[i] + ep_n[i] + dg * G[i];
  }
  for (int i=0; i<nh; i++) {
    R[i+6] = -alpha_np1[i] + alpha_n[i] + dg * G[i + 6];
  }

  // Check convergence
  double Rn = norm2_vec(R, nhist());
  if (verbose_) {
    std::cout << "\t" << Rn << "\t" << f << std::endl;
  }
  if ((Rn < tol_) && (f < tol_)) return true;

  // Call down for jacobian
  double J[nhist()*nhist()];
  get_jacobian(s, q, Di, T, mu, dg, J);

  // We need to do a few solves, the first two can't be in-place
  int n = nhist();
  double Rsolve[nhist()];
  std::copy(R, R+nhist(), Rsolve);
  solve_mat(J, n, Rsolve);
  
  double Gsolve[nhist()];
  std::copy(G, G+nhist(), Gsolve);
  solve_mat(J, n, Gsolve);

  // Calculate the increment
  double ddg = (f - dot_vec(G, Rsolve, nhist()))/dot_vec(G, Gsolve, nhist());

  // Calculate the increments
  for (int i=0; i<nhist(); i++) R[i] += ddg * G[i];
  solve_mat(J, n, R);

  for (int i=0; i<6; i++) {
    ep_np1[i] += R[i] / (2.0 * mu);
  }
  for (int i=0; i<nh; i++) {
    for (int j=0; j<nh; j++) {
      alpha_np1[i] += Di[CINDEX(i,j,nh)] * R[i + 6]; 
    }
  }

  dg += ddg;

  return false;
}

void RIAFModel::get_jacobian(const double * const s, const double * const q,
                             const double * const Di,
                             double T, double mu, double dg, double * const J) const
{
  std::fill(J, J+nhist()*nhist(), 0.0);

  // No choice but to block I think
  size_t nh = hardening_->nhist();
  double J11[6*6];
  surface_->df_dsds(s, q, T, J11);

  for (int i=0; i<6*6; i++) J11[i] *= dg;
  for (int i=0; i<6; i++) J11[CINDEX(i,i,6)] += 1.0 / (2.0 * mu);

  double J12[6*nh];
  surface_->df_dsdq(s, q, T, J12);
  for (int i=0; i<6*nh; i++) J12[i] *= dg;

  double J21[nh*6];
  surface_->df_dqds(s, q, T, J21);
  for (int i=0; i<nh*6; i++) J21[i] *= dg;

  double J22[nh*nh];
  surface_->df_dqdq(s, q, T, J22);
  for (int i=0; i<nh*nh; i++) J22[i] = J22[i] * dg + Di[i];

  // Stick everything in place
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      J[CINDEX(i,j,nhist())] = J11[CINDEX(i,j,6)];
    }
  }
  for (int i=0; i<6; i++) {
    for (int j=0; j<nh; j++) {
      J[CINDEX(i,j+6,nhist())] = J12[CINDEX(i,j,nh)];
      J[CINDEX(j+6,i,nhist())] = J21[CINDEX(j,i,6)];
    }
  }
  for (int i=0; i<nh; i++) {
    for (int j=0; j<nh; j++) {
      J[CINDEX(i+6,j+6,nhist())] = J22[CINDEX(i,j,nh)];
    }
  }
  
  /*
  std::cout << "CHECKING" << std::endl;
  for (int i=0; i<nhist(); i++) {
    for (int j=0; j<nhist(); j++) {
      std::cout << J[CINDEX(i,j,nhist())] << "\t";
    }
    std::cout << std::endl;
  }
  */
}

void RIAFModel::get_tangent(double mu, double T, 
                            const double * const alpha_np1, double dg,
                            const double * const s_np1, 
                            double * const A_np1) const
{
  double mod[36];  
  if (dg > 0.0) {
    // Get to stress space
    int nh = hardening_->nhist();
    double q_np1[nh];
    hardening_->q(alpha_np1, T, q_np1);

    // Grab required derivatives
    double grad[6];
    
    surface_->df_ds(s_np1, q_np1, T, grad);
    surface_->df_dsds(s_np1, q_np1, T, mod);

    for (int i=0; i<36; i++) {
      mod[i] *= dg;
    }
    for (int i=0; i<6; i++) {
      mod[CINDEX(i,i,6)] += 1.0 / (2.0 * mu);
    }

    // Invert in place
    invert_mat(mod, 6);

    // Setup outer
    double N[6];
    mat_vec(mod, 6, grad, 6, N);
    double v = std::sqrt(dot_vec(N, grad, 6));
    for (int i=0; i<6; i++) {
      N[i] /= v;
    }

    outer_update_minus(N, 6, N, 6, mod);
  }
  else {
    std::fill(mod,mod+36,0.0);
    for (int i=0; i<6; i++) {
      mod[CINDEX(i,i,6)] = 2.0 * mu;
    }
  }

  // Now take care of the deviatoric part
  double cr[36] = {
    2.0/3.0, -1.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
    -1.0/3.0, 2.0/3.0, -1.0/3.0, 0.0, 0.0, 0.0,
    -1.0/3.0, -1.0/3.0, 2.0/3.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

  mat_mat(6, 6, 6, mod, cr, A_np1);

}


} // namespace neml
