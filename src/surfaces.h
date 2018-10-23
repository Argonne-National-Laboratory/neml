#ifndef SURFACES_H
#define SURFACES_H

#include <cstddef>
#include <stdarg.h>
#include <memory>
#include <algorithm>
#include <string>

#include "objects.h"
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
class YieldSurface: public NEMLObject {
 public:
  /// Indicates how many history variables the model expects to get
  virtual size_t nhist() const = 0;

  /// Yield function
  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const = 0;

  /// Gradient wrt stress
  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const = 0;
  /// Gradient wrt history
  virtual int df_dq(const double* const s, const double* const q, double T,
                double * const df) const = 0;

  /// Hessian dsds
  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const = 0;
  /// Hessian dqdq
  virtual int df_dqdq(const double* const s, const double* const q, double T,
                double * const ddf) const = 0;
  /// Hessian dsdq
  virtual int df_dsdq(const double* const s, const double* const q, double T,
                double * const ddf) const = 0;
  /// Hessian dqds
  virtual int df_dqds(const double* const s, const double* const q, double T,
                double * const ddf) const = 0;
};

/// Helper to reduce a isotropic + kinematic function to isotropic only
template<class BT, typename... Args>
class IsoFunction: public YieldSurface {
 public:
  /// Take whatever args the template class takes
  IsoFunction(Args... args) :
      base_(new BT(args...))
  {

  }
  
  /// Also interfaces with a single isotropic hardening variable
  virtual size_t nhist() const
  {
    return 1;
  }
  
  /// Just call with zero kinematic hardening
  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const
  {
    double * qn = expand_hist_(q);
    int ier = base_->f(s, qn, T, fv);
    delete [] qn;
    return ier;
  }
  
  /// Call with zero kinematic hardening
  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const
  {
    double * qn = expand_hist_(q);
    int ier = base_->df_ds(s, qn, T, df);
    delete [] qn;
    return ier;
  }

  /// Call with zero kinematic hardening
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

  /// Call with zero kinematic hardening
  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const
  {
    double * qn = expand_hist_(q);
    int ier = base_->df_dsds(s, qn, T, ddf);
    delete [] qn;
    return ier;
  }

  /// Call with zero kinematic hardening
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

  /// Call with zero kinematic hardening
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

  /// Call with zero kinematic hardening
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
  /// No parameters
  IsoKinJ2();
  virtual ~IsoKinJ2();
  
  /// String type for object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();
 
  /// Expects 7 history variables [isotropic 6-Mandel-vector-backstress]
  virtual size_t nhist() const;
  
  /// J2(stress + backstress) + sqrt(2/3) * isotropic
  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const;
  
  /// Gradient wrt stress
  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const;
  /// Gradient wrt history
  virtual int df_dq(const double* const s, const double* const q, double T,
                double * const df) const;
  
  /// Hessian dsds
  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const;
  /// Hessian dqdq
  virtual int df_dqdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  /// Hessian dsdq
  virtual int df_dsdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  /// Hessian dqds
  virtual int df_dqds(const double* const s, const double* const q, double T,
                double * const ddf) const;
 
};

static Register<IsoKinJ2> regIsoKinJ2;

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
  /// No parameters
  IsoJ2() :
      IsoFunction<IsoKinJ2>()
  {
  }
  
  /// String type for object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();
};

static Register<IsoJ2> regIsoJ2;

/// Combined isotropic/kinematic hardening with some mean stress contribution
//
//  History variables are:
//    hist[0]     q (isotropic hardening)
//    hist[1:7]   X (backstress)
//
class IsoKinJ2I1: public YieldSurface {
 public:
  /// Parameters: h prefactor, l exponent
  IsoKinJ2I1(std::shared_ptr<Interpolate> h, std::shared_ptr<Interpolate> l);
  virtual ~IsoKinJ2I1();
  
  /// String type for object system
  static std::string type();
  /// Initialize from parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Expects 7 history variables [isotropic 6-Mandel-vector-backstress]
  virtual size_t nhist() const;
  
  /// J2(stress + backstress) + isotropic + sign(mean_stress) * h *
  /// |mean_stress|^l
  virtual int f(const double* const s, const double* const q, double T,
                double & fv) const;
  
  /// Gradient wrt stress
  virtual int df_ds(const double* const s, const double* const q, double T,
                double * const df) const;
  /// Gradient wrt q
  virtual int df_dq(const double* const s, const double* const q, double T,
                double * const df) const;
  
  /// Hessian dsds
  virtual int df_dsds(const double* const s, const double* const q, double T,
                double * const ddf) const;
  /// Hessian dqdq
  virtual int df_dqdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  /// Hessian dsdq
  virtual int df_dsdq(const double* const s, const double* const q, double T,
                double * const ddf) const;
  /// Hessian dqds
  virtual int df_dqds(const double* const s, const double* const q, double T,
                double * const ddf) const;

 private:
  const std::shared_ptr<Interpolate> h_;
  const std::shared_ptr<Interpolate> l_;
 
};

static Register<IsoKinJ2I1> regIsoKinJ2I1;

/// Isotropic only version of J2I1 surface
class IsoJ2I1: public IsoFunction<IsoKinJ2I1, std::shared_ptr<Interpolate>,
    std::shared_ptr<Interpolate>> {
 public:
  // h prefactor and l exponent
  IsoJ2I1(std::shared_ptr<Interpolate> h, std::shared_ptr<Interpolate> l) :
      IsoFunction<IsoKinJ2I1, std::shared_ptr<Interpolate>, 
      std::shared_ptr<Interpolate>>(h, l)
  {
  }

  /// String type for object system
  static std::string type();
  /// Initialize from parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

};

static Register<IsoJ2I1> regIsoJ2I1;

} // namespace neml

#endif // SURFACES_H
