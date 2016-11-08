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

} // namespace neml

#endif // SURFACES_H
