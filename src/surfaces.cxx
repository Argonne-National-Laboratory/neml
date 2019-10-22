#include "surfaces.h"

#include "math/nemlmath.h"
#include "nemlerror.h"

#include <limits>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace neml {

std::string IsoJ2::type()
{
  return "IsoJ2";
}

std::unique_ptr<NEMLObject> IsoJ2::initialize(ParameterSet & params)
{
  return neml::make_unique<IsoJ2>(); 
}

ParameterSet IsoJ2::parameters()
{
  ParameterSet pset(IsoJ2::type());

  return pset;
}

// Combined isotropic/kinematic J2 surface
IsoKinJ2::IsoKinJ2()
{

}

IsoKinJ2::~IsoKinJ2()
{

}

std::string IsoKinJ2::type()
{
  return "IsoKinJ2";
}

std::unique_ptr<NEMLObject> IsoKinJ2::initialize(ParameterSet & params)
{
  return neml::make_unique<IsoKinJ2>(); 
}

ParameterSet IsoKinJ2::parameters()
{
  ParameterSet pset(IsoKinJ2::type());

  return pset;
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

  if (nv > 0.0) {
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

  if (nv > 0.0) {
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
  
  if (nv > 0.0) {
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
  }

  return 0;
}

int IsoKinJ2::df_dqds(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  std::fill(ddf, ddf+nhist()*6, 0.0);
  
  double ss[36];
  int ier = df_dsds(s, q, T, ss);
  if (ier != SUCCESS) return ier;

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      ddf[CINDEX((i+1),j,6)] = ss[CINDEX(i,j,6)];
    }
  }

  return 0;
}

// J2I1 implementation
// Combined isotropic/kinematic with some I1 contribution
IsoKinJ2I1::IsoKinJ2I1(std::shared_ptr<Interpolate> h, 
                       std::shared_ptr<Interpolate> l) :
    h_(h), l_(l)
{

}

IsoKinJ2I1::~IsoKinJ2I1()
{

}

std::string IsoKinJ2I1::type()
{
  return "IsoKinJ2I1";
}

std::unique_ptr<NEMLObject> IsoKinJ2I1::initialize(ParameterSet & params)
{
  return neml::make_unique<IsoKinJ2I1>(
      params.get_object_parameter<Interpolate>("h"),
      params.get_object_parameter<Interpolate>("l")
      ); 
}

ParameterSet IsoKinJ2I1::parameters()
{
  ParameterSet pset(IsoKinJ2I1::type());

  pset.add_parameter<NEMLObject>("h");
  pset.add_parameter<NEMLObject>("l");

  return pset;
}

size_t IsoKinJ2I1::nhist() const
{
  return 7;
}

int IsoKinJ2I1::f(const double* const s, const double* const q, double T,
              double & fv) const
{
  double sdev[6];
  std::copy(s, s+6, sdev);
  dev_vec(sdev);
  add_vec(sdev, &q[1], 6, sdev);
  fv = norm2_vec(sdev, 6) + sqrt(2.0/3.0) * q[0] + 
      copysign(h_->value(T) * pow( fabs(s[0] + s[1] + s[2]), l_->value(T) ),
               s[0] + s[1] + s[2]);
  return 0;
}

int IsoKinJ2I1::df_ds(const double* const s, const double* const q, double T,
              double * const df) const
{
  std::copy(s, s+6, df);
  dev_vec(df);
  add_vec(df, &q[1], 6, df);
  normalize_vec(df, 6);

  // Compute dsh/ds
  double dsh[6];
  for (int i=0; i<3; i++) {
    dsh[i] = h_->value(T) * l_->value(T) * pow( fabs(s[0] + s[1] + s[2]),
                                               l_->value(T) - 1.0 );
  }
  for (int i=3; i<6; i++) {
    dsh[i] = 0.0;
  }
  add_vec(df, dsh, 6, df);

  return 0;
}

int IsoKinJ2I1::df_dq(const double* const s, const double* const q, double T,
              double * const df) const
{
  df[0] = sqrt(2.0/3.0);
  std::copy(s, s+6, &df[1]);
  dev_vec(&df[1]);
  add_vec(&df[1], &q[1], 6, &df[1]);
  normalize_vec(&df[1], 6);

  return 0;
}

int IsoKinJ2I1::df_dsds(const double* const s, const double* const q, double T,
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
  
  // Compute ddsh/dsds
  double iv2[6];
  for (int i=0; i<3; i++) {
    iv2[i] = copysign(h_->value(T) * l_->value(T) * (l_->value(T) - 1.0) *
                      pow( fabs(s[0]+s[1]+s[2]), l_->value(T) - 2.0),
                      s[0] + s[1] + s[2]);    
  }
  for (int i=3; i<6; i++) {
    iv2[i] = 0.0;
  }
  outer_update(iv2, 6, jv, 6, ddf);

  return 0;
}

int IsoKinJ2I1::df_dqdq(const double* const s, const double* const q, double T,
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

int IsoKinJ2I1::df_dsdq(const double* const s, const double* const q, double T,
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

int IsoKinJ2I1::df_dqds(const double* const s, const double* const q, double T,
              double * const ddf) const
{
  std::fill(ddf, ddf+nhist()*6, 0.0);
  
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

  outer_update_minus(iv, 6, jv, 6, ss);

  outer_update_minus(n, 6, n, 6, ss);
  for (int i=0; i<36; i++) {
    ss[i] /= nv;
  }

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      ddf[CINDEX((i+1),j,6)] = ss[CINDEX(i,j,6)];
    }
  }

  return 0;
}

std::string IsoJ2I1::type()
{
  return "IsoJ2I1";
}

std::unique_ptr<NEMLObject> IsoJ2I1::initialize(ParameterSet & params)
{
  return neml::make_unique<IsoJ2I1>(
      params.get_object_parameter<Interpolate>("h"),
      params.get_object_parameter<Interpolate>("l")
      ); 
}

ParameterSet IsoJ2I1::parameters()
{
  ParameterSet pset(IsoJ2I1::type());

  pset.add_parameter<NEMLObject>("h");
  pset.add_parameter<NEMLObject>("l");

  return pset;
}


} // namespace neml
