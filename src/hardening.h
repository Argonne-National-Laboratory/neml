#ifndef HARDENING_H
#define HARDENING_H

#include <cstddef>

namespace neml {

/// Interface for a generic hardening rule
//    1) Take alpha to q
//    2) Give the gradient of that function
class HardeningRule {
 public:
  virtual size_t nhist() const = 0;
  virtual int init_hist(double * const alpha) const = 0;
  virtual int q(const double * const alpha, double T, double * const qv) const = 0;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const = 0;
};

/// Isotropic hardening rules
class IsotropicHardeningRule: public HardeningRule {
 public:
  virtual size_t nhist() const;
  virtual int init_hist(double * const alpha) const;
  virtual int q(const double * const alpha, double T, double * const qv) const = 0;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const = 0;
};

/// Linear, isotropic hardening
class LinearIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  LinearIsotropicHardeningRule(double s0, double K);
  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

  double s0() const;
  double K() const;
 private:
  const double s0_, K_;
};

} // namespace neml

#endif // HARDENING_H
