#include "hardening.h"

#include "nemlmath.h"

namespace neml {

AssociativeHardening::AssociativeHardening()
{

}

AssociativeHardening::~AssociativeHardening()
{

}

int AssociativeHardening::D_inv(const double * const alpha, double T, double * const Dv) const
{
  D(alpha, T, Dv);
  return invert_mat(Dv, nhist());
}

IsoJ2LinearAHardening::IsoJ2LinearAHardening(double K0, double Kp) :
    AssociativeHardening(), K0_(K0), Kp_(Kp)
{

}

IsoJ2LinearAHardening::~IsoJ2LinearAHardening()
{

}

size_t IsoJ2LinearAHardening::nhist() const
{
  return 1;
}

int IsoJ2LinearAHardening::init_hist(double * const alpha) const
{
  alpha[0] = 0.0;
  return 0;
}

int IsoJ2LinearAHardening::q(const double * const alpha, double T, double * const qv) const
{
  qv[0] = -K0_ - Kp_ * alpha[0];
  return 0;
}

int IsoJ2LinearAHardening::D(const double * const alpha, double T, double * const Dv) const
{
  Dv[0] = -Kp_;
  return 0;
}

int IsoJ2LinearAHardening::D_inv(const double * const alpha, double T, double * const Dv) const
{
  Dv[0] = -1.0 / Kp_;
}

double IsoJ2LinearAHardening::K0() const
{
  return K0_;
}


double IsoJ2LinearAHardening::Kp() const
{
  return Kp_;
}

} // namespace neml
