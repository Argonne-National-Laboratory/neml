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

} // namespace neml
