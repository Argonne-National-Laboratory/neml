#include "surfaces.h"
#include "nemlmath.h"

#include <limits>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace neml {

// Pure isotropic J2 surface
IsoJ2::IsoJ2()
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
  double sdev[6];
  std::copy(s, s+6, sdev);
  dev_vec(sdev);
  fv = norm2_vec(sdev, 6) + sqrt(2.0/3.0) * q[0];
  return 0;
}

int IsoJ2::df_ds(const double* const s, const double* const q, double T,
              double * const df) const
{
  std::copy(s, s+6, df);
  dev_vec(df);
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
  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  double nv = norm2_vec(n, 6);
  normalize_vec(n, 6);
  
  std::fill(ddf, ddf+36, 0.0);
  for (int i=0; i<6; i++) {
    ddf[CINDEX(i,i,6)] += 1.0;
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

  outer_update_minus(iv, 6, jv, 6, ddf);

  outer_update_minus(n, 6, n, 6, ddf);
  for (int i=0; i<36; i++) {
    ddf[i] /= nv;
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


// Combined isotropic/kinematic J2 surface
IsoKinJ2::IsoKinJ2()
{

}

IsoKinJ2::~IsoKinJ2()
{

}

size_t IsoKinJ2::nhist() const
{
  return 7;
}

int IsoKinJ2::f(const double* const s, const double* const q, double T,
              double & fv) const
{
  double sdev[6];
  std::copy(s, s+6, sdev);
  dev_vec(sdev);
  add_vec(sdev, &q[1], 6, sdev);
  fv = norm2_vec(sdev, 6) + sqrt(2.0/3.0) * q[0];
  return 0;
}

int IsoKinJ2::df_ds(const double* const s, const double* const q, double T,
              double * const df) const
{
  std::copy(s, s+6, df);
  dev_vec(df);
  add_vec(df, &q[1], 6, df);
  normalize_vec(df, 6);
  return 0;
}

int IsoKinJ2::df_dq(const double* const s, const double* const q, double T,
              double * const df) const
{
  df[0] = sqrt(2.0/3.0);
  df_ds(s, q, T, &df[1]);

  return 0;
}

int IsoKinJ2::df_dsds(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, &q[1], 6, n);
  double nv = norm2_vec(n, 6);
  normalize_vec(n, 6);
  
  std::fill(ddf, ddf+36, 0.0);
  for (int i=0; i<6; i++) {
    ddf[CINDEX(i,i,6)] += 1.0;
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

  outer_update_minus(iv, 6, jv, 6, ddf);

  outer_update_minus(n, 6, n, 6, ddf);
  for (int i=0; i<36; i++) {
    ddf[i] /= nv;
  }

  return 0;
}

int IsoKinJ2::df_dqdq(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  std::fill(ddf, ddf+nhist()*nhist(), 0.0);
  
  double ss[36];
  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, &q[1], 6, n);
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

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      ddf[CINDEX((i+1),(j+1),nhist())] = ss[CINDEX(i,j,6)];
    }
  }

  return 0;
}

int IsoKinJ2::df_dsdq(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  std::fill(ddf, ddf+6*nhist(), 0.0);
  
  double ss[36];
  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, &q[1], 6, n);
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

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      ddf[CINDEX(i,(j+1),nhist())] = ss[CINDEX(i,j,6)];
    }
  }

  return 0;
}

int IsoKinJ2::df_dqds(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  std::fill(ddf, ddf+nhist()*6, 0.0);
  
  double ss[36];
  df_dsds(s, q, T, ss);

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      ddf[CINDEX((i+1),j,6)] = ss[CINDEX(i,j,6)];
    }
  }

  return 0;
}


} // namespace neml
