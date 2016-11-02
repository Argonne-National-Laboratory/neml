#ifndef SURFACE_H
#define SURFACE_H

#include <cstddef>

namespace neml {

class YieldSurface {
 public:
  YieldSurface();
  virtual ~YieldSurface();
  
  // Interface
  virtual size_t nhist() const = 0;

  // Yield function
  virtual int f(const double* const s, const double* const h, double T,
                double & fv) const = 0;

  // Gradients
  virtual int df_ds(const double* const s, const double* const h, double T,
                double * const df) const = 0;
  virtual int df_dq(const double* const s, const double* const h, double T,
                double * const df) const = 0;

  // Hessian
  virtual int df_dsds(const double* const s, const double* const h, double T,
                double * const ddf) const = 0;
  virtual int df_dqdq(const double* const s, const double* const h, double T,
                double * const ddf) const = 0;
  virtual int df_dsdq(const double* const s, const double* const h, double T,
                double * const ddf) const = 0;
  virtual int df_dqds(const double* const s, const double* const h, double T,
                double * const ddf) const = 0;
};

/// Just isotropic hardening with a von Mises surface
//
//  History variables are:
//    hist[0]     K (isotropic hardening)
//
//
class IsoJ2: public YieldSurface {
 public:
  IsoJ2();
  virtual ~IsoJ2();
 
  // Defined interface
  virtual size_t nhist() const;

  virtual int f(const double* const s, const double* const h, double T,
                double & fv) const;

  virtual int df_ds(const double* const s, const double* const h, double T,
                double * const df) const;
  virtual int df_dq(const double* const s, const double* const h, double T,
                double * const df) const;

  virtual int df_dsds(const double* const s, const double* const h, double T,
                double * const ddf) const;
  virtual int df_dqdq(const double* const s, const double* const h, double T,
                double * const ddf) const;
  virtual int df_dsdq(const double* const s, const double* const h, double T,
                double * const ddf) const;
  virtual int df_dqds(const double* const s, const double* const h, double T,
                double * const ddf) const;
 
};

/// Kinematic and isotropic hardening with the usual von Mises surface
//
//  History variables are:
//    hist[0]     K (isotropic hardening parameter)
//    hist[1:7]   X (backstress)
//
//
class KinIsoJ2: public YieldSurface {
 public:
  KinIsoJ2();
  virtual ~KinIsoJ2();
 
  // Defined interface
  virtual size_t nhist() const;

  virtual int f(const double* const s, const double* const h, double T,
                double & fv) const;

  virtual int df_ds(const double* const s, const double* const h, double T,
                double * const df) const;
  virtual int df_dq(const double* const s, const double* const h, double T,
                double * const df) const;

  virtual int df_dsds(const double* const s, const double* const h, double T,
                double * const ddf) const;
  virtual int df_dqdq(const double* const s, const double* const h, double T,
                double * const ddf) const;
  virtual int df_dsdq(const double* const s, const double* const h, double T,
                double * const ddf) const;
  virtual int df_dqds(const double* const s, const double* const h, double T,
                double * const ddf) const;
 
};

} // namespace neml

#endif // SURFACE_H
