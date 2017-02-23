#include "hardening.h"

#include "nemlmath.h"

#include <cmath>
#include <algorithm>

namespace neml {



size_t IsotropicHardeningRule::nhist() const
{
  return 1;
}

int IsotropicHardeningRule::init_hist(double * const alpha) const
{
  alpha[0] = 0.0;

  return 0;
}

// Implementation of linear hardening
LinearIsotropicHardeningRule::LinearIsotropicHardeningRule(double s0, double K) :
    s0_(new ConstantInterpolate(s0)), K_(new ConstantInterpolate(K))
{

}

LinearIsotropicHardeningRule::LinearIsotropicHardeningRule(std::shared_ptr<Interpolate> s0, std::shared_ptr<Interpolate> K) :
    s0_(s0), K_(K)
{

}

int LinearIsotropicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  qv[0] = -s0_->value(T) - K_->value(T) * alpha[0];

  return 0;
}

int LinearIsotropicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  dqv[0] = -K_->value(T);

  return 0;
}

double LinearIsotropicHardeningRule::s0(double T) const
{
  return s0_->value(T);
}

double LinearIsotropicHardeningRule::K(double T) const
{
  return K_->value(T);
}

// Implementation of voce hardening
VoceIsotropicHardeningRule::VoceIsotropicHardeningRule(double s0, double R,
                                                       double d) :
    s0_(new ConstantInterpolate(s0)), R_(new ConstantInterpolate(R)),
    d_(new ConstantInterpolate(d))
{

}

VoceIsotropicHardeningRule::VoceIsotropicHardeningRule(std::shared_ptr<Interpolate> s0, std::shared_ptr<Interpolate> R, std::shared_ptr<Interpolate> d) :
    s0_(s0), R_(R), d_(d)
{

}


int VoceIsotropicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  qv[0] = -s0_->value(T) - R_->value(T) * (1.0 - exp(-d_->value(T) * alpha[0]));

  return 0;
}

int VoceIsotropicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  dqv[0] = -d_->value(T) * R_->value(T) * exp(-d_->value(T) * alpha[0]);

  return 0;
}

double VoceIsotropicHardeningRule::s0(double T) const
{
  return s0_->value(T);
}

double VoceIsotropicHardeningRule::R(double T) const
{
  return R_->value(T);
}

double VoceIsotropicHardeningRule::d(double T) const
{
  return d_->value(T);
}

// Implementation of kinematic base class
size_t KinematicHardeningRule::nhist() const
{
  return 6;
}

int KinematicHardeningRule::init_hist(double * const alpha) const
{
  for (int i=0; i<6; i++) alpha[0] = 0.0;

  return 0;
}

// Implementation of linear kinematic hardening
LinearKinematicHardeningRule::LinearKinematicHardeningRule(double H) :
    H_(new ConstantInterpolate(H))
{

}

LinearKinematicHardeningRule::LinearKinematicHardeningRule(std::shared_ptr<Interpolate> H) :
    H_(H)
{

}

int LinearKinematicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  for (int i=0; i<6; i++) {
    qv[i] = -H_->value(T) * alpha[i];
  }

  return 0;
}

int LinearKinematicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  std::fill(dqv, dqv+36, 0.0);
  for (int i=0; i<6; i++) {
    dqv[CINDEX(i,i,6)] = -H_->value(T);
  }

  return 0;
}

double LinearKinematicHardeningRule::H(double T) const
{
  return H_->value(T);
}

CombinedHardeningRule::CombinedHardeningRule(
    std::shared_ptr<IsotropicHardeningRule> iso,
    std::shared_ptr<KinematicHardeningRule> kin) :
      iso_(iso), kin_(kin)
{

}
size_t CombinedHardeningRule::nhist() const
{
  return iso_->nhist() + kin_->nhist();
}

int CombinedHardeningRule::init_hist(double * const alpha) const
{
  int ier = iso_->init_hist(alpha);
  return kin_->init_hist(&alpha[iso_->nhist()]);
}

int CombinedHardeningRule::q(const double * const alpha, double T, 
                             double * const qv) const
{
  iso_->q(alpha, T, qv);
  return kin_->q(&alpha[iso_->nhist()], T, &qv[iso_->nhist()]);
}

int CombinedHardeningRule::dq_da(const double * const alpha, double T, 
                                 double * const dqv) const
{
  // Annoying this doesn't work nicely...
  double id[iso_->nhist()*iso_->nhist()];
  int ier = iso_->dq_da(alpha, T, id);
  double kd[kin_->nhist()*kin_->nhist()];
  ier = kin_->dq_da(&alpha[iso_->nhist()], T, kd);

  std::fill(dqv, dqv + nhist()*nhist(), 0.0);

  for (int i=0; i<iso_->nhist(); i++) {
    for (int j=0; j<iso_->nhist(); j++) {
      dqv[CINDEX(i,j,nhist())] = id[CINDEX(i,j,iso_->nhist())];
    }
  }
  
  int os = iso_->nhist();
  for (int i=0; i<kin_->nhist(); i++) {
    for (int j=0; j<kin_->nhist(); j++) {
      dqv[CINDEX((i+os),(j+os),nhist())] = kd[CINDEX(i,j,kin_->nhist())];
    }
  }

  return 0;
}


// Provide zeros for these
int NonAssociativeHardening::h_time(const double * const s, 
                                    const double * const alpha, double T,
                                    double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
  return 0;
}

int NonAssociativeHardening::dh_ds_time(const double * const s,
                                        const double * const alpha, double T,
                                        double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*6, 0.0);
  return 0;
}

int NonAssociativeHardening::dh_da_time(const double * const s,
                                        const double * const alpha, double T,
                                        double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);
  return 0;
}

// Provide zeros for these
int NonAssociativeHardening::h_temp(const double * const s,
                                    const double * const alpha, double T,
                                    double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
  return 0;
}

int NonAssociativeHardening::dh_ds_temp(const double * const s,
                                        const double * const alpha, double T,
                                        double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*6, 0.0);
  return 0;
}

int NonAssociativeHardening::dh_da_temp(const double * const s,
                                        const double * const alpha, double T,
                                        double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);
  return 0;
}

// Begin non-associative hardening rules
//
// Gamma functions for Chaboche
//
// Constant
//
ConstantGamma::ConstantGamma(double g) :
    g_(new ConstantInterpolate(g))
{

}

ConstantGamma::ConstantGamma(std::shared_ptr<Interpolate> g) :
    g_(g)
{

}

double ConstantGamma::gamma(double ep, double T) const {
  return g_->value(T);
}

double ConstantGamma::dgamma(double ep, double T) const {
  return 0;
}

double ConstantGamma::g(double T) const {
  return g_->value(T);
}

//
// Saturating
//
SatGamma::SatGamma(double gs, double g0, double beta) :
    gs_(new ConstantInterpolate(gs)), g0_(new ConstantInterpolate(g0)),
    beta_(new ConstantInterpolate(beta))
{

}

SatGamma::SatGamma(std::shared_ptr<Interpolate> gs, std::shared_ptr<Interpolate> g0, std::shared_ptr<Interpolate> beta) :
    gs_(gs), g0_(g0), beta_(beta)
{

}

double SatGamma::gamma(double ep, double T) const {
  return gs_->value(T) + (g0_->value(T) - gs_->value(T)) * exp(-beta_->value(T) * ep);
}

double SatGamma::dgamma(double ep, double T) const {
  return beta_->value(T) * (gs_->value(T) - g0_->value(T)) * exp(-beta_->value(T) * ep);
}

double SatGamma::gs(double T) const {
  return gs_->value(T);
}

double SatGamma::g0(double T) const {
  return g0_->value(T);
}

double SatGamma::beta(double T) const {
  return beta_->value(T);
}

//
// Chaboche
//
Chaboche::Chaboche(std::shared_ptr<IsotropicHardeningRule> iso,
                   std::vector<double> c,
                   std::vector<std::shared_ptr<GammaModel>> gmodels) :
    iso_(iso), c_(make_vector(c)), gmodels_(gmodels), n_(c.size()),
    relax_(false)
{

}

Chaboche::Chaboche(std::shared_ptr<IsotropicHardeningRule> iso,
           std::vector<std::shared_ptr<Interpolate>> c,
           std::vector<std::shared_ptr<GammaModel>> gmodels) :
    iso_(iso), c_(c), gmodels_(gmodels), n_(c.size()), relax_(false)
{

}

Chaboche::Chaboche(std::shared_ptr<IsotropicHardeningRule> iso,
                   std::vector<double> c,
                   std::vector<std::shared_ptr<GammaModel>> gmodels,
                   std::vector<double> A,
                   std::vector<double> a) :
    iso_(iso), c_(make_vector(c)), gmodels_(gmodels), n_(c.size()),
    A_(make_vector(A)), a_(make_vector(a)), relax_(true)
{

}

Chaboche::Chaboche(std::shared_ptr<IsotropicHardeningRule> iso,
           std::vector<std::shared_ptr<Interpolate>> c,
           std::vector<std::shared_ptr<GammaModel>> gmodels,
           std::vector<std::shared_ptr<Interpolate>> A,
           std::vector<std::shared_ptr<Interpolate>> a) :
    iso_(iso), c_(c), gmodels_(gmodels), n_(c.size()), A_(A), a_(a), 
    relax_(true)
{

}

size_t Chaboche::ninter() const
{
  return 1 + 6;
}

size_t Chaboche::nhist() const
{
  return 1 + 6 * n_;
}

int Chaboche::init_hist(double * const alpha) const
{
  std::fill(alpha, alpha+nhist(), 0.0);
  return 0;
}

int Chaboche::q(const double * const alpha, double T, double * const qv) const
{
  iso_->q(alpha, T, qv);
  std::fill(qv+1, qv+7, 0.0);
  
  // Helps with unrolling
  int n = n_;

  for (int i=0; i<n; i++) {
    for (int j=0; j<6; j++) {
      qv[j+1] += alpha[1+i*6+j];
    }
  }
  return 0;
}

int Chaboche::dq_da(const double * const alpha, double T, double * const qv) const
{
  std::fill(qv, qv+(ninter()*nhist()), 0.0);
  iso_->dq_da(alpha, T, qv); // fills in (0,0)
  
  // Help unroll
  int n = n_;

  for (int i=0; i<n; i++) {
    for (int j=0; j<6; j++) {
      qv[CINDEX((j+1),(1+i*6+j), nhist())] = 1.0;
    }
  }
  return 0;
}

int Chaboche::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  hv[0] = sqrt(2.0/3.0); // Isotropic part

  double X[6], nv[6];
  backstress_(alpha, X);
  std::copy(s, s+6, nv);
  dev_vec(nv);
  add_vec(nv, X, 6, nv);
  normalize_vec(nv, 6);
  
  // Note the extra factor of sqrt(2.0/3.0) -- this is to make it equivalent
  // to Chaboche's original definition
  
  std::vector<double> c = eval_vector(c_, T);

  for (int i=0; i<n_; i++) {
    for (int j=0; j<6; j++) {
      hv[1+i*6+j] = - 2.0 / 3.0 * c[i] * nv[j] - sqrt(2.0/3.0) * gmodels_[i]->gamma(alpha[0], T) * alpha[1+i*6+j];
    }
  }

  return 0;
}

int Chaboche::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv + nhist()*6, 0.0);

  std::vector<double> c = eval_vector(c_, T);

  double X[6];
  backstress_(alpha, X);

  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, X, 6, n);
  double nv = norm2_vec(n, 6);
  normalize_vec(n, 6);
  
  double nn[36];

  std::fill(nn, nn+36, 0.0);
  for (int i=0; i<6; i++) {
    nn[CINDEX(i,i,6)] += 1.0;
  }
  
  double iv[6];
  double jv[6];
  for (int i=0; i<3; i++) {
    iv[i] = 1.0 / 3.0;
    jv[i] = 1.0;
  }
  for (int i=3; i<6; i++) {
    iv[i] = 0.0;
    jv[i] = 0.0;
  }

  outer_update_minus(iv, 6, jv, 6, nn);

  outer_update_minus(n, 6, n, 6, nn);
  for (int i=0; i<36; i++) {
    nn[i] /= nv;
  }
  
  // Fill in...
  for (int i=0; i<n_; i++) {
    for (int j=0; j<6; j++) {
      for (int k=0; k<6; k++) {
        dhv[CINDEX((1+i*6+j),(k),6)] = -2.0 / 3.0 * c[i] * nn[CINDEX(j,k,6)];
      }
    }
  }
  
  return 0;
}

int Chaboche::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv + nhist()*nhist(), 0.0);

  std::vector<double> c = eval_vector(c_, T);

  double X[6];
  backstress_(alpha, X);

  double ss[36];
  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, X, 6, n);
  double nv = norm2_vec(n, 6);
  normalize_vec(n, 6);
  
  std::fill(ss, ss+36, 0.0);
  for (int i=0; i<6; i++) {
    ss[CINDEX(i,i,6)] += 1.0;
  }
  
  outer_update_minus(n, 6, n, 6, ss);
  for (int i=0; i<36; i++) {
    ss[i] /= nv;
  }
  
  // Fill in the gamma part
  for (int i=0; i<n_; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((1+i*6+j),(1+i*6+j),nhist())] -= sqrt(2.0/3.0) * gmodels_[i]->gamma(alpha[0], T);
    }
  }

  // Fill in the ss part
  for (int bi=0; bi<n_; bi++) {
    for (int i=0; i<6; i++) {
      for (int bj=0; bj<n_; bj++) {
        for (int j=0; j<6; j++) {
          dhv[CINDEX((1+bi*6+i),(1+bj*6+j),nhist())] -= 2.0 / 3.0 * c[bi]  * 
              ss[CINDEX(i,j,6)];
        }
      }
    }
  }

  // Fill in the alpha part
  for (int i=0; i<n_; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((1+i*6+j),0,nhist())] = -sqrt(2.0/3.0) * 
          gmodels_[i]->dgamma(alpha[0], T) * alpha[1+i*6+j];
    }
  }

  return 0;
}

int Chaboche::h_time(const double * const s, const double * const alpha, 
                     double T, double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
  if (not relax_) return 0;
 
  std::vector<double> A = eval_vector(A_, T);
  std::vector<double> a = eval_vector(a_, T);

  double Xi[6];
  double nXi;
  for (int i=0; i<n_; i++) {
    std::copy(&alpha[1+i*6], &alpha[1+(i+1)*6], Xi);
    nXi = norm2_vec(Xi, 6);
    for (int j=0; j<6; j++) {
      hv[1+i*6+j] = -A[i] * sqrt(3.0/2.0) * pow(nXi, a[i] - 1.0) *
          alpha[1+i*6+j];
    }
  }

  return 0;
}

int Chaboche::dh_ds_time(const double * const s, const double * const alpha, 
                         double T, double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*6, 0.0);
  if (not relax_) return 0;

  // Also return if relax

  return 0;
}

int Chaboche::dh_da_time(const double * const s, const double * const alpha,
                         double T, double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);
  if (not relax_) return 0;

  std::vector<double> A = eval_vector(A_, T);
  std::vector<double> a = eval_vector(a_, T);

  int nh = nhist();
  int n = n_;
  
  double XX[36];
  double Xi[6];
  double nXi;
  int ia,ib;
  double d;
  for (int i=0; i<n; i++) {
    std::copy(&alpha[1+i*6], &alpha[1+(i+1)*6], Xi);
    nXi = norm2_vec(Xi, 6);
    normalize_vec(Xi, 6);
    outer_vec(Xi, 6, Xi, 6, XX);
    for (int j=0; j<6; j++) {
      ia = 1 + i*6 + j;
      for (int k=0; k<6; k++) {
        ib = 1 + i*6 + k;
        if (j == k) {
          d = 1.0;
        }
        else {
          d = 0.0;
        }
        dhv[CINDEX(ia,ib,nh)] = -A[i] * sqrt(3.0/2.0) * pow(nXi, a[i]-1.0) * (
            d + (a[i] - 1.0) * XX[CINDEX(j,k,6)]);
      }
    }
  }

  return 0;
}

int Chaboche::n() const
{
  return n_;
}

std::vector<double> Chaboche::c(double T) const
{
  return eval_vector(c_, T);
}

void Chaboche::backstress_(const double * const alpha, double * const X) const
{
  std::fill(X, X+6, 0.0);
  for (int i=0; i<n_; i++) {
    for (int j=0; j<6; j++) {
      X[j] += alpha[1+i*6+j];
    }
  }
}

} // namespace neml
