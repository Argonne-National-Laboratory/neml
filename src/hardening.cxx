#include "hardening.h"

#include <cmath>

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

} // namespace neml
