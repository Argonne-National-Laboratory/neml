#ifndef SURFACES_H
#define SURFACES_H

#include <cstddef>

namespace neml {

class YieldSurface {
 public:
  // Interface
  virtual size_t nhist() const = 0;

  // Yield function
  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const = 0;

  // Gradients
  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const = 0;
  virtual int df_dq(const double* const s, const double* const q, double T,
                double * const df) const = 0;

  // Hessian
  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const = 0;
  virtual int df_dqdq(const double* const s, const double* const q, double T,
                double * const ddf) const = 0;
  virtual int df_dsdq(const double* const s, const double* const q, double T,
                double * const ddf) const = 0;
  virtual int df_dqds(const double* const s, const double* const q, double T,
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

  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const;

  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const;
  virtual int df_dq(const double* const s, const double* const q, double T,
                double * const df) const;

  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dqdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dsdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dqds(const double* const s, const double* const q, double T,
                double * const ddf) const;
 
};

/// Combined isotropic/kinematic hardening with a von Mises surface
//
//  History variables are:
//    hist[0]     q (isotropic hardening)
//    hist[1:7]   X (backstress)
//
class IsoKinJ2: public YieldSurface {
 public:
  IsoKinJ2();
  virtual ~IsoKinJ2();
 
  // Defined interface
  virtual size_t nhist() const;

  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const;

  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const;
  virtual int df_dq(const double* const s, const double* const q, double T,
                double * const df) const;

  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dqdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dsdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dqds(const double* const s, const double* const q, double T,
                double * const ddf) const;
 
};

/// Combined isotropic/kinematic hardening with some mean stress contribution
//
//  History variables are:
//    hist[0]     q (isotropic hardening)
//    hist[1:7]   X (backstress)
//
class IsoKinJ2I1: public YieldSurface {
 public:
  IsoKinJ2I1(const double h, const double l);
  virtual ~IsoKinJ2I1();
 
  // Defined interface
  virtual size_t nhist() const;

  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const;

  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const;
  virtual int df_dq(const double* const s, const double* const q, double T,
                double * const df) const;

  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dqdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dsdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  virtual int df_dqds(const double* const s, const double* const q, double T,
                double * const ddf) const;

 private:
  const double h_, l_;
 
};

} // namespace neml

#endif // SURFACES_H
