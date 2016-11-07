#include "surface.h"
#include "nemlmath.h"

#include <limits>
#include <algorithm>

namespace neml {

// Base class implementation

YieldSurface::YieldSurface()
{

}

YieldSurface::~YieldSurface()
{

}

// Pure isotropic J2 surface
IsoJ2::IsoJ2() :
    YieldSurface()
{

}

IsoJ2::~IsoJ2()
{

}

size_t IsoJ2::nhist() const
{
  return 1;
}

int IsoJ2::f(const double* const s, const double* const q, double T,
              double & fv) const
{
  fv = norm2_vec(s, 6) + sqrt(2.0/3.0) * q[0];
  return 0;
}

int IsoJ2::df_ds(const double* const s, const double* const q, double T,
              double * const df) const
{
  std::copy(s, s+6, df);
  normalize_vec(df, 6);

  return 0;
}

int IsoJ2::df_dq(const double* const s, const double* const q, double T,
              double * const df) const
{
  df[0] = sqrt(2.0/3.0);

  return 0;
}

int IsoJ2::df_dsds(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  double s1[6], s2[6];
  std::copy(s, s+6, s1);

  double ns = norm2_vec(s1, 6);

  if (ns < std::numeric_limits<double>::epsilon()) {
    std::fill(ddf, ddf+36, 0.0);
    return 0;
  }

  normalize_vec(s1, 6);
  std::copy(s1, s1+6, s2);
  minus_vec(s2, 6);
  for (int i=0; i<6; i++) {
    s2[i] /= ns;
  }

  outer_vec(s1, 6, s2, 6, ddf);

  for (int i=0; i<6; i++) {
    ddf[CINDEX(i,i,6)] += 1.0 / ns;
  }

  return 0;
}

int IsoJ2::df_dqdq(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  std::fill(ddf, ddf+nhist()*nhist(), 0.0);

  return 0;
}

int IsoJ2::df_dsdq(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  std::fill(ddf, ddf+6*nhist(), 0.0);

  return 0;
}

int IsoJ2::df_dqds(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  std::fill(ddf, ddf+nhist()*6, 0.0);

  return 0;
}



// Kinematic and isotropic hardening J2 surface
KinIsoJ2::KinIsoJ2() :
    YieldSurface()
{

}

KinIsoJ2::~KinIsoJ2()
{

}

size_t KinIsoJ2::nhist() const
{
  return 7;
}

int KinIsoJ2::f(const double* const s, const double* const q, double T,
              double & fv) const
{

}

int KinIsoJ2::df_ds(const double* const s, const double* const q, double T,
              double * const df) const
{

}

int KinIsoJ2::df_dq(const double* const s, const double* const q, double T,
              double * const df) const
{

}

int KinIsoJ2::df_dsds(const double* const s, const double* const q, double T,
              double * const ddf) const
{

}

int KinIsoJ2::df_dqdq(const double* const s, const double* const q, double T,
              double * const ddf) const
{

}

int KinIsoJ2::df_dsdq(const double* const s, const double* const q, double T,
              double * const ddf) const
{

}

int KinIsoJ2::df_dqds(const double* const s, const double* const q, double T,
              double * const ddf) const
{

}


} // namespace neml
