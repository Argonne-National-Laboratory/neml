#include "hardening.h"

namespace neml {

size_t IsotropicHardeningRule::nhist() const
{
  return 1;
}

int IsotropicHardeningRule::init_hist(double * const alpha) const
{
  alpha[0] = 0.0;
}

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
}

double LinearIsotropicHardeningRule::s0() const
{
  return s0_;
}

double LinearIsotropicHardeningRule::K() const
{
  return K_;
}

} // namespace neml
