#include "creep.h"

#include "nemlmath.h"

#include <cmath>

namespace neml {

// Scalar creep default derivatives for time and temperature
int ScalarCreepRule::dg_dt(double seq, double eeq, double t, double T,
                           double & dg) const
{
  dg = 0.0;
  return 0;
}

int ScalarCreepRule::dg_dT(double seq, double eeq, double t, double T,
                           double & dg) const
{
  dg = 0.0;
  return 0;
}

PowerLawCreep::PowerLawCreep(double A, double n) :
    A_(new ConstantInterpolate(A)), n_(new ConstantInterpolate(n))
{

}

PowerLawCreep::PowerLawCreep(std::shared_ptr<Interpolate> A,
                             std::shared_ptr<Interpolate> n) :
    A_(A), n_(n)
{

}

int PowerLawCreep::g(double seq, double eeq, double t, double T, double & g) const
{
  g = A_->value(T) * pow(seq, n_->value(T));
  return 0;
}

int PowerLawCreep::dg_ds(double seq, double eeq, double t, double T, double & dg) const
{
  double nv = n_->value(T);

  dg = A_->value(T) * nv * pow(seq, nv - 1.0);
  return 0;
}

int PowerLawCreep::dg_de(double seq, double eeq, double t, double T, double & dg) const
{
  dg = 0.0;
  return 0;
}

double PowerLawCreep::A(double T) const
{
  return A_->value(T);
}

double PowerLawCreep::n(double T) const
{
  return n_->value(T);
}

} // namespace neml
