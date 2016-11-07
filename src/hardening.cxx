#include "hardening.h"

#include "nemlmath.h"
#include <cmath>

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
  Dv[0] = Kp_;
  return 0;
}

int IsoJ2LinearAHardening::D_inv(const double * const alpha, double T, double * const Dv) const
{
  Dv[0] = 1.0 / Kp_;
}

double IsoJ2LinearAHardening::K0() const
{
  return K0_;
}


double IsoJ2LinearAHardening::Kp() const
{
  return Kp_;
}


IsoJ2VoceAHardening::IsoJ2VoceAHardening(double K0, double Ksat, double delta) :
    AssociativeHardening(), K0_(K0), Ksat_(Ksat), delta_(delta)
{

}

IsoJ2VoceAHardening::~IsoJ2VoceAHardening()
{

}

size_t IsoJ2VoceAHardening::nhist() const
{
  return 1;
}

int IsoJ2VoceAHardening::init_hist(double * const alpha) const
{
  alpha[0] = 0.0;
  return 0;
}

int IsoJ2VoceAHardening::q(const double * const alpha, double T, double * const qv) const
{
  qv[0] = -K0_ - Ksat_ * (1.0 - std::exp(-delta_ * alpha[0]));
  return 0;
}

int IsoJ2VoceAHardening::D(const double * const alpha, double T, double * const Dv) const
{
  Dv[0] = delta_ * Ksat_ * std::exp(-delta_ * alpha[0]);
  return 0;
}

int IsoJ2VoceAHardening::D_inv(const double * const alpha, double T, double * const Dv) const
{
  Dv[0] = 1.0 / (delta_ * Ksat_ * std::exp(-delta_ * alpha[0]));
}

double IsoJ2VoceAHardening::K0() const
{
  return K0_;
}

double IsoJ2VoceAHardening::Ksat() const
{
  return Ksat_;
}

double IsoJ2VoceAHardening::delta() const
{
  return delta_;
}

} // namespace neml
