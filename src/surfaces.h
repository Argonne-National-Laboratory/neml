#ifndef SURFACES_H
#define SURFACES_H

#include <cstddef>
#include <stdarg.h>
#include <memory>
#include <algorithm>

#include "nemlmath.h"
#include "interpolate.h"

namespace neml {

/// Base class for generic yield surfaces
//  Yield surfaces are defined in terms of the yield function, 1st
//  derivatives wrt. the stress and history, and second derivatives
//  wrt. the same.
//
//  Hardening models must match the number of history variables but
//  right now there's no enforcement mechanism to make sure that a model
//  with the correct number of variables also has the right type of
//  variables.
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

/// Helper to reduce a isotropic + kinematic function to isotropic only
template<class BT, typename... Args>
class IsoFunction: public YieldSurface {
 public:
  IsoFunction(Args... args) :
      base_(new BT(args...))
  {

  }

  virtual size_t nhist() const
  {
    return 1;
  }

  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const
  {
    double * qn = expand_hist_(q);
    int ier = base_->f(s, qn, T, fv);
    delete [] qn;
    return ier;
  }

  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const
  {
    double * qn = expand_hist_(q);
    int ier = base_->df_ds(s, qn, T, df);
    delete [] qn;
    return ier;
  }

  virtual int df_dq(const double* const s, const double* const q, double T,
                double * const df) const
  {
    double * qn = expand_hist_(q);
    double * dfn = new double[base_->nhist()];
    int ier = base_->df_dq(s, qn, T, dfn);
    df[0] = dfn[0];
    delete [] qn;
    delete [] dfn;
    return ier;
  }

  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const
  {
    double * qn = expand_hist_(q);
    int ier = base_->df_dsds(s, qn, T, ddf);
    delete [] qn;
    return ier;
  }

  virtual int df_dqdq(const double* const s, const double* const q, double T,
                double * const ddf) const
  {
    double * qn = expand_hist_(q);
    double * ddfn = new double[(base_->nhist())*(base_->nhist())];
    int ier = base_->df_dqdq(s, qn, T, ddfn);
    ddf[0] = ddfn[0];
    delete [] qn;
    delete [] ddfn;
    return ier;
  }

  virtual int df_dsdq(const double* const s, const double* const q, double T,
                double * const ddf) const
  {
    // This one is annoying
    double * qn = expand_hist_(q);
    double * ddfn = new double[6*(base_->nhist())];
    int ier = base_->df_dsdq(s, qn, T, ddfn);
    for (int i=0; i<6; i++) {
      ddf[i] = ddfn[CINDEX(i,0,base_->nhist())];
    }
    delete [] qn;
    delete [] ddfn;
    return ier;
  }

  virtual int df_dqds(const double* const s, const double* const q, double T,
                double * const ddf) const
  {
    double * qn = expand_hist_(q);
    double * ddfn = new double[(base_->nhist())*6];
    int ier = base_->df_dqds(s, qn, T, ddfn);
    std::copy(ddfn,ddfn+6,ddf);
    delete [] qn;
    delete [] ddfn;
    return ier;
  }

 private:
  double * expand_hist_(const double* const q) const 
  {
    double * qn = new double[7];
    qn[0] = q[0];
    std::fill(qn+1,qn+7,0.0);
    return qn;
  }

 private:
  std::unique_ptr<BT> base_;

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


/// Just isotropic hardening with a von Mises surface
//
//  History variables are:
//    hist[0]     K (isotropic hardening)
//
//  Originally I coded a custom version of this that didn't re-use the
//  IsoKinJ2 code.  I switched to this for convenience and reliability but
//  note it is slightly memory inefficient.
//
class IsoJ2: public IsoFunction<IsoKinJ2> {
 public:
  IsoJ2() :
      IsoFunction<IsoKinJ2>()
  {
  }
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
  IsoKinJ2I1(std::shared_ptr<Interpolate> h, std::shared_ptr<Interpolate> l);
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
  const std::shared_ptr<Interpolate> h_;
  const std::shared_ptr<Interpolate> l_;
 
};

/// Isotropic only version of J2I1 surface
class IsoJ2I1: public IsoFunction<IsoKinJ2I1, std::shared_ptr<Interpolate>,
    std::shared_ptr<Interpolate>> {
 public:
  IsoJ2I1(std::shared_ptr<Interpolate> h, std::shared_ptr<Interpolate> l) :
      IsoFunction<IsoKinJ2I1, std::shared_ptr<Interpolate>, 
      std::shared_ptr<Interpolate>>(h, l)
  {
  }

};

} // namespace neml

#endif // SURFACES_H
