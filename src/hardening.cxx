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
    s0_(s0), K_(K)
{

}

int LinearIsotropicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  qv[0] = -s0_ - K_ * alpha[0];

  return 0;
}

int LinearIsotropicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  dqv[0] = -K_;

  return 0;
}

double LinearIsotropicHardeningRule::s0() const
{
  return s0_;
}

double LinearIsotropicHardeningRule::K() const
{
  return K_;
}

// Implementation of voce hardening
VoceIsotropicHardeningRule::VoceIsotropicHardeningRule(double s0, double R,
                                                       double d) :
    s0_(s0), R_(R), d_(d)
{

}

int VoceIsotropicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  qv[0] = -s0_ - R_ * (1.0 - exp(-d_ * alpha[0]));

  return 0;
}

int VoceIsotropicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  dqv[0] = -d_ * R_ * exp(-d_ * alpha[0]);

  return 0;
}

double VoceIsotropicHardeningRule::s0() const
{
  return s0_;
}

double VoceIsotropicHardeningRule::R() const
{
  return R_;
}

double VoceIsotropicHardeningRule::d() const
{
  return d_;
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
    H_(H)
{

}

int LinearKinematicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  for (int i=0; i<6; i++) {
    qv[i] = -H_ * alpha[i];
  }

  return 0;
}

int LinearKinematicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  std::fill(dqv, dqv+36, 0.0);
  for (int i=0; i<6; i++) {
    dqv[CINDEX(i,i,6)] = -H_;
  }

  return 0;
}

double LinearKinematicHardeningRule::H() const
{
  return H_;
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



// Begin non-associative hardening rules
//
// Chaboche
//
Chaboche::Chaboche(std::shared_ptr<IsotropicHardeningRule> iso,
                   int n, const double * const c, const double * const r)
  : iso_(iso), n_(n), c_(c, c+n), r_(r, r+n)
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
  for (int i=0; i<n_; i++) {
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
  for (int i=0; i<n_; i++) {
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

  for (int i=0; i<n_; i++) {
    for (int j=0; j<6; j++) {
      hv[1+i*6+j] = - c_[i] * r_[i] * (nv[j] - alpha[1+i*6+j] / r_[i]);
    }
  }

  return 0;
}

int Chaboche::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv + nhist()*6, 0.0);

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
        dhv[CINDEX((1+i*6+j),(k),6)] = -c_[i] * r_[i] * nn[CINDEX(j,k,6)];
      }
    }
  }
  
  return 0;
}

int Chaboche::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv + nhist()*nhist(), 0.0);

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
  
  // Fill in the ci part
  for (int i=0; i<n_; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((1+i*6+j),(1+i*6+j),nhist())] += c_[i];
    }
  }

  // Fill in the ss part
  for (int bi=0; bi<n_; bi++) {
    for (int i=0; i<6; i++) {
      for (int bj=0; bj<n_; bj++) {
        for (int j=0; j<6; j++) {
          dhv[CINDEX((1+bi*6+i),(1+bj*6+j),nhist())] -= c_[bi] * r_[bi] * 
              ss[CINDEX(i,j,6)];
        }
      }
    }
  }

  return 0;
}

int Chaboche::n() const
{
  return n_;
}

const std::vector<double> & Chaboche::c() const
{
  return c_;
}

const std::vector<double> & Chaboche::r() const
{
  return r_;
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
