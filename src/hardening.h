#ifndef HARDENING_H
#define HARDENING_H

#include "interpolate.h"

#include <cstddef>
#include <memory>
#include <vector>

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

  double s0(double T) const;
  double K(double T) const;
 private:
  const std::shared_ptr<const Interpolate> s0_, K_;
};

/// Voce isotropic hardening
class VoceIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  VoceIsotropicHardeningRule(double s0, double R, double d);
  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

  double s0(double T) const;
  double R(double T) const;
  double d(double T) const;

 private:
  const std::shared_ptr<const Interpolate> s0_, R_, d_;

};

class KinematicHardeningRule: public HardeningRule {
 public:
  virtual size_t nhist() const;
  virtual int init_hist(double * const alpha) const;
  virtual int q(const double * const alpha, double T, double * const qv) const = 0;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const = 0;
};

class LinearKinematicHardeningRule: public KinematicHardeningRule {
 public:
  LinearKinematicHardeningRule(double H);
  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

  double H(double T) const;

 private:
  const std::shared_ptr<const Interpolate> H_;
};

class CombinedHardeningRule: public HardeningRule {
 public:
  CombinedHardeningRule(std::shared_ptr<IsotropicHardeningRule> iso,
                        std::shared_ptr<KinematicHardeningRule> kin);
  virtual size_t nhist() const;
  virtual int init_hist(double * const alpha) const;
  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

 private:
  std::shared_ptr<IsotropicHardeningRule> iso_;
  std::shared_ptr<KinematicHardeningRule> kin_;
};

/// ABC of a non-associative hardening rule
class NonAssociativeHardening {
 public:
  virtual size_t ninter() const = 0; // How many "q" variables it spits out
  virtual size_t nhist() const = 0; // How many internal variables it stores

  virtual int init_hist(double * const alpha) const = 0;

  virtual int q(const double * const alpha, double T, double * const qv) const = 0;
  virtual int dq_da(const double * const alpha, double T, double * const qv) const = 0;

  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const = 0;
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
};

/// Gamma models for the Chaboche backstress
class GammaModel {
 public:
  virtual double gamma(double ep, double T) const = 0;
  virtual double dgamma(double ep, double T) const = 0;

};

class ConstantGamma: public GammaModel {
 public:
  ConstantGamma(double g);

  virtual double gamma(double ep, double T) const;
  virtual double dgamma(double ep, double T) const;

  double g(double T) const;

 private:
  const std::shared_ptr<const Interpolate> g_;

};

class SatGamma: public GammaModel {
 public:
  SatGamma(double gs, double g0, double beta);

  virtual double gamma(double ep, double T) const;
  virtual double dgamma(double ep, double T) const;

  double gs(double T) const;
  double g0(double T) const;
  double beta(double T) const;

 private:
  const std::shared_ptr<const Interpolate> gs_, g0_, beta_;

};

/// Chaboche model: generalized Frederick-Armstrong
//    This model degenerates to Frederick-Armstrong for n = 1
class Chaboche: public NonAssociativeHardening {
 public:
  /// New interface with vectors
  Chaboche(std::shared_ptr<IsotropicHardeningRule> iso,
           std::vector<double> c,
           std::vector<std::shared_ptr<GammaModel>> gmodels);

  /// Older interface assuming constant gammas
  Chaboche(std::shared_ptr<IsotropicHardeningRule> iso,
                     int n, const double * const c, const double * const r);

  virtual size_t ninter() const; // How many "q" variables it spits out
  virtual size_t nhist() const; // How many internal variables it stores

  virtual int init_hist(double * const alpha) const;

  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const qv) const;

  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  // Getters
  int n() const;
  std::vector<double> c(double T) const;

 private:
  void backstress_(const double * const alpha, double * const X) const;

  std::shared_ptr<IsotropicHardeningRule> iso_;
  const int n_;
  const std::vector<std::shared_ptr<const Interpolate>> c_;
  std::vector<std::shared_ptr<GammaModel>> gmodels_;
};


} // namespace neml

#endif // HARDENING_H
