#include "surface.h"
#include "nemlmath.h"

namespace neml {

// Base class implementation

YieldSurface::YieldSurface()
{

}

YieldSurface::~YieldSurface()
{

}

int YieldSurface::D_inv(const double* const s, const double* const h, double T,
                double * const Dv) const
{
  D(s, h, T, Dv);
  return invert_matrix(Dv, nhist());
}


// Kinematic and isotropic hardening ala Simo and Hughes
KinIsoJ2::KinIsoJ2() :
    YieldSurface()
{

}

KinIsoJ2::~KinIsoJ2()
{

}

size_t KinIsoJ2::nhist() const
{
  return 2;
}

int KinIsoJ2::f(const double* const s, const double* const h, double T,
              double & fv) const
{

}

int KinIsoJ2::df_ds(const double* const s, const double* const h, double T,
              double * const df) const
{

}

int KinIsoJ2::df_dq(const double* const s, const double* const h, double T,
              double * const df) const
{

}

int KinIsoJ2::df_dsds(const double* const s, const double* const h, double T,
              double * const ddf) const
{

}

int KinIsoJ2::df_dqdq(const double* const s, const double* const h, double T,
              double * const ddf) const
{

}

int KinIsoJ2::df_dsdq(const double* const s, const double* const h, double T,
              double * const ddf) const
{

}

int KinIsoJ2::df_dqds(const double* const s, const double* const h, double T,
              double * const ddf) const
{

}

int KinIsoJ2::D(const double* const s, const double* const h, double T,
              double * const Dv) const
{

}

int KinIsoJ2::D_inv(const double* const s, const double* const h, double T,
              double * const Dv) const
{

}

// Linear hardening model
LinearKinIsoJ2::LinearKinIsoJ2(double K0, double Kb, double Hb) :
    KinIsoJ2(), K0_(K0), Kb_(Kb), Hb_(Hb)
{

}

LinearKinIsoJ2::~LinearKinIsoJ2()
{

}

double LinearKinIsoJ2::K(double a) const
{
  return K0_ + Kb_ * a;
}

double LinearKinIsoJ2::Kp(double a) const
{
  return Kb_;
}

double LinearKinIsoJ2::Hp(double a) const
{
  return Hb_;
}

// Getters
double LinearKinIsoJ2::K0() const
{
  return K0_;
}

double LinearKinIsoJ2::Kb() const
{
  return Kb_;
}

double LinearKinIsoJ2::Hb() const
{
  return Hb_;
}


} // namespace neml
