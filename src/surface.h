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

  // Generalized plastic moduli
  virtual int D(const double* const s, const double* const h, double T,
                double * const Dv) const = 0;
  // Actually more common thing to need, so give the option to override
  // By default invert with lapack
  virtual int D_inv(const double* const s, const double* const h, double T,
                double * const Dv) const;

};

/// Kinematic and isotropic hardening with the usual von Mises surface
//
//  Hardening described by functions K(a) and H'(a)
//
class KinIsoJ2: public YieldSurface {
 public:
  KinIsoJ2();
  virtual ~KinIsoJ2();
 
  // New abstract interface
  virtual double K(double a) const = 0;
  virtual double Kp(double a) const = 0;
  virtual double Hp(double a) const = 0;

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

  virtual int D(const double* const s, const double* const h, double T,
                double * const Dv) const;
  virtual int D_inv(const double* const s, const double* const h, double T,
                double * const Dv) const;
};

// A simple implementation of the isotropic and kinematic hardening:
// K = K0 + Kb * a, K' = Kb
// Hp = Hb
class LinearKinIsoJ2: public KinIsoJ2 {
 public:
  LinearKinIsoJ2(double K0, double Kb, double Hb);
  virtual ~LinearKinIsoJ2();
  
  // Implementation
  virtual double K(double a) const;
  virtual double Kp(double a) const;
  virtual double Hp(double a) const;
  
  double K0() const;
  double Kb() const;
  double Hb() const;

 private:
  const double K0_, Kb_, Hb_;

};


} // namespace neml

#endif // SURFACE_H
