#ifndef HARDENING_H
#define HARDENING_H

#include "interpolate.h"
#include "objects.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace neml {

/// Interface for a generic hardening rule
//    1) Take alpha to q
//    2) Give the gradient of that function
class HardeningRule: public NEMLObject {
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
  LinearIsotropicHardeningRule(std::shared_ptr<Interpolate> s0, std::shared_ptr<Interpolate> K);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

  double s0(double T) const;
  double K(double T) const;
 private:
  const std::shared_ptr<const Interpolate> s0_, K_;
};

static Register<LinearIsotropicHardeningRule> regLinearIsotropicHardeningRule;

/// Isotropic hardening with flow stress from some interpolation function
//    The convention will be to provide a flow curve as
//    (plastic strain, flow stress) tuples, with the value of the curve at
//    0 being the initial yield stress
class InterpolatedIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  InterpolatedIsotropicHardeningRule(std::shared_ptr<Interpolate> flow);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

 private:
  const std::shared_ptr<const Interpolate> flow_;

};

static Register<InterpolatedIsotropicHardeningRule>
  regInterpolatedIsotropicHardeningRule;

/// Voce isotropic hardening
class VoceIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  VoceIsotropicHardeningRule(std::shared_ptr<Interpolate> s0, std::shared_ptr<Interpolate> R, std::shared_ptr<Interpolate> d);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

  double s0(double T) const;
  double R(double T) const;
  double d(double T) const;

 private:
  const std::shared_ptr<const Interpolate> s0_, R_, d_;
};

static Register<VoceIsotropicHardeningRule> regVoceIsotropicHardeningRule;

/// Combined hardening rule superimposing a bunch of separate ones
class CombinedIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  CombinedIsotropicHardeningRule(
      std::vector<std::shared_ptr<IsotropicHardeningRule>> rules);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

  size_t nrules() const;

 private:
  const std::vector<std::shared_ptr<IsotropicHardeningRule>> rules_;

};

static Register<CombinedIsotropicHardeningRule> regCombinedIsotropicHardeningRule;

/// Base class for pure kinematic hardening
class KinematicHardeningRule: public HardeningRule {
 public:
  virtual size_t nhist() const;
  virtual int init_hist(double * const alpha) const;
  virtual int q(const double * const alpha, double T, double * const qv) const = 0;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const = 0;
};

/// Simple linear kinematic hardening
class LinearKinematicHardeningRule: public KinematicHardeningRule {
 public:
  LinearKinematicHardeningRule(std::shared_ptr<Interpolate> H);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

  double H(double T) const;

 private:
  const std::shared_ptr<const Interpolate> H_;
};

static Register<LinearKinematicHardeningRule> regLinearKinematicHardeningRule;

/// Class to combine isotropic and kinematic hardening rules
class CombinedHardeningRule: public HardeningRule {
 public:
  CombinedHardeningRule(std::shared_ptr<IsotropicHardeningRule> iso,
                        std::shared_ptr<KinematicHardeningRule> kin);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual size_t nhist() const;
  virtual int init_hist(double * const alpha) const;
  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int dq_da(const double * const alpha, double T, double * const dqv) const;

 private:
  std::shared_ptr<IsotropicHardeningRule> iso_;
  std::shared_ptr<KinematicHardeningRule> kin_;
};

static Register<CombinedHardeningRule> regCombinedHardeningRule;

/// ABC of a non-associative hardening rule
class NonAssociativeHardening: public NEMLObject {
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

  // Hardening rule wrt to time
  virtual int h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  // Hardening rule wrt to temperature
  virtual int h_temp(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
};

/// Gamma as a function of equivalent plastic strain
class GammaModel: public NEMLObject {
 public:
  virtual double gamma(double ep, double T) const = 0;
  virtual double dgamma(double ep, double T) const = 0;

};

/// Gamma is just a consant
class ConstantGamma: public GammaModel {
 public:
  ConstantGamma(std::shared_ptr<Interpolate> g);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual double gamma(double ep, double T) const;
  virtual double dgamma(double ep, double T) const;

  double g(double T) const;

 private:
  const std::shared_ptr<const Interpolate> g_;
};

static Register<ConstantGamma> regConstantGamma;

/// Gamma evolves with a saturating Voce form
class SatGamma: public GammaModel {
 public:
  SatGamma(std::shared_ptr<Interpolate> gs, std::shared_ptr<Interpolate> g0, std::shared_ptr<Interpolate> beta);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

  virtual double gamma(double ep, double T) const;
  virtual double dgamma(double ep, double T) const;

  double gs(double T) const;
  double g0(double T) const;
  double beta(double T) const;

 private:
  const std::shared_ptr<const Interpolate> gs_, g0_, beta_;
};

static Register<SatGamma> regSatGamma;

/// Chaboche model: generalized Frederick-Armstrong
//    This model degenerates to Frederick-Armstrong for n = 1
class Chaboche: public NonAssociativeHardening {
 public:
  Chaboche(std::shared_ptr<IsotropicHardeningRule> iso,
           std::vector<std::shared_ptr<Interpolate>> c,
           std::vector<std::shared_ptr<GammaModel>> gmodels,
           std::vector<std::shared_ptr<Interpolate>> A,
           std::vector<std::shared_ptr<Interpolate>> a,
           bool noniso = true);

  static std::string type();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();

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

  // Hardening rule wrt to time
  virtual int h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  
  // Hardening rule wrt temperature
  virtual int h_temp(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  // Getters
  int n() const;
  std::vector<double> c(double T) const;

 private:
  void backstress_(const double * const alpha, double * const X) const;
  
 private:
  std::shared_ptr<IsotropicHardeningRule> iso_;
  const int n_;
  const std::vector<std::shared_ptr<Interpolate>> c_;
  std::vector<std::shared_ptr<GammaModel>> gmodels_;

  const std::vector<std::shared_ptr<Interpolate>> A_;
  const std::vector<std::shared_ptr<Interpolate>> a_;
  const bool relax_;

  const bool noniso_;
};

static Register<Chaboche> regChaboche;

} // namespace neml

#endif // HARDENING_H
