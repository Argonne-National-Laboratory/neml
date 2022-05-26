#ifndef HARDENING_H
#define HARDENING_H

#include "interpolate.h"
#include "objects.h"
#include "history.h"

#include "windows.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace neml {

/// Interface for a generic hardening rule
//    1) Take alpha to q
//    2) Give the gradient of that function
class NEML_EXPORT HardeningRule: public HistoryNEMLObject {
 public:
  HardeningRule(ParameterSet & params);

  /// The map between strain-like and stress-like variables
  virtual void q(const double * const alpha, double T, double * const qv) const = 0;
  /// The derivative of the map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const = 0;
};

/// Isotropic hardening rules
class NEML_EXPORT IsotropicHardeningRule: public HardeningRule {
 public:
  IsotropicHardeningRule(ParameterSet & params);
  /// Setup the internal state
  virtual void populate_hist(History & h) const;
  /// Initialize the history
  virtual void init_hist(History & h) const;
  /// Map between the strain variables and the single isotropic hardening
  /// parameter
  virtual void q(const double * const alpha, double T, double * const qv) const = 0;
  /// Derivative of the map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const = 0;
};

/// Linear, isotropic hardening
class NEML_EXPORT LinearIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  /// Parameters: initial surface size and linear coefficient
  LinearIsotropicHardeningRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// q = -s0 - K * alpha[0]
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const;

  /// Getter for the yield stress
  double s0(double T) const;

  /// Getter for the linear coefficient
  double K(double T) const;
 private:
  const std::shared_ptr<const Interpolate> s0_, K_;
};

static Register<LinearIsotropicHardeningRule> regLinearIsotropicHardeningRule;

/// Isotropic hardening with flow stress from some interpolation function
//    The convention will be to provide a flow curve as
//    (plastic strain, flow stress) tuples, with the value of the curve at
//    0 being the initial yield stress
class NEML_EXPORT InterpolatedIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  /// Parameter is the interpolate to use
  InterpolatedIsotropicHardeningRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// q = -interpolate(T)
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of the map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const;

 private:
  const std::shared_ptr<const Interpolate> flow_;
};

static Register<InterpolatedIsotropicHardeningRule>
  regInterpolatedIsotropicHardeningRule;

/// Voce isotropic hardening
class NEML_EXPORT VoceIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  /// Parameters: initial yield stress, total increase amount,
  /// saturation speed constant
  VoceIsotropicHardeningRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// q = -s0 - R * (1 - exp(-d * alpha[0]))
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const;

  /// Getter for initial yield stress
  double s0(double T) const;
  /// Getter for R
  double R(double T) const;
  /// Getter for d
  double d(double T) const;

 private:
  const std::shared_ptr<const Interpolate> s0_, R_, d_;
};

static Register<VoceIsotropicHardeningRule> regVoceIsotropicHardeningRule;

/// Power law hardening
class NEML_EXPORT PowerLawIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  /// Parameters: initial yield stress, prefactor, exponent
  PowerLawIsotropicHardeningRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// q = -s0 - A * alpha[0]**n
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const;

 private:
  const std::shared_ptr<const Interpolate> s0_, A_, n_;
};

static Register<PowerLawIsotropicHardeningRule> regPowerLawIsotropicHardeningRule;

/// Combined hardening rule superimposing a bunch of separate ones
class NEML_EXPORT CombinedIsotropicHardeningRule: public IsotropicHardeningRule {
 public:
  /// Parameter is a vector of isotropic hardening rules
  CombinedIsotropicHardeningRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// q = Sum(q_i(alpha))
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const;

  /// Getter on the number of combined rules
  size_t nrules() const;

 private:
  const std::vector<std::shared_ptr<IsotropicHardeningRule>> rules_;

};

static Register<CombinedIsotropicHardeningRule> regCombinedIsotropicHardeningRule;

/// Base class for pure kinematic hardening
class NEML_EXPORT KinematicHardeningRule: public HardeningRule {
 public:
  KinematicHardeningRule(ParameterSet & params);
  /// Setup the internal state
  virtual void populate_hist(History & h) const;
  /// Initialize the history
  virtual void init_hist(History & h) const;

  /// Map between the backstrain and the backstress
  virtual void q(const double * const alpha, double T, double * const qv) const = 0;
  /// Derivative of the map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const = 0;
};

/// Simple linear kinematic hardening
class NEML_EXPORT LinearKinematicHardeningRule: public KinematicHardeningRule {
 public:
  /// Parameter is the linear hardening coefficient
  LinearKinematicHardeningRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// q = -H * alpha[:6]
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of the map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const;

  /// Getter for the hardening coefficieint
  double H(double T) const;

 private:
  const std::shared_ptr<const Interpolate> H_;
};

static Register<LinearKinematicHardeningRule> regLinearKinematicHardeningRule;

/// Class to combine isotropic and kinematic hardening rules
class NEML_EXPORT CombinedHardeningRule: public HardeningRule {
 public:
  /// Parameters: a isotropic hardening rule and a kinematic hardening rule
  CombinedHardeningRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup the internal state
  virtual void populate_hist(History & h) const;
  /// Initialize the history
  virtual void init_hist(History & h) const;
  /// q[0] = isotropic model, q[1:7] = kinematic model
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of the map
  virtual void dq_da(const double * const alpha, double T, double * const dqv) const;

 private:
  std::shared_ptr<IsotropicHardeningRule> iso_;
  std::shared_ptr<KinematicHardeningRule> kin_;
};

static Register<CombinedHardeningRule> regCombinedHardeningRule;

/// ABC of a non-associative hardening rule
class NEML_EXPORT NonAssociativeHardening: public HistoryNEMLObject {
 public:
  NonAssociativeHardening(ParameterSet & params);
  /// How many stress-like variables
  virtual size_t ninter() const = 0; // How many "q" variables it spits out
  
  /// Map from strain to stress
  virtual void q(const double * const alpha, double T, double * const qv) const = 0;
  /// Derivative of the map
  virtual void dq_da(const double * const alpha, double T, double * const qv) const = 0;

  /// Hardening proportional to the equivalent inelastic strain
  virtual void h(const double * const s, const double * const alpha, double T,
                double * const hv) const = 0;
  /// Derivative of h wrt stress
  virtual void dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
  /// Derivative of h wrt history
  virtual void dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;

  /// Hardening proportional to time
  virtual void h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_time wrt stress
  virtual void dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_time wrt history
  virtual void dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening proportional to the temperature rate
  virtual void h_temp(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_temp wrt stress
  virtual void dh_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_temp wrt history
  virtual void dh_da_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
};

/// Model for the gamma constant used in Chaboche hardening
class NEML_EXPORT GammaModel: public NEMLObject {
 public:
  GammaModel(ParameterSet & params);
  /// Gamma as a function of equivalent inelastic strain
  virtual double gamma(double ep, double T) const = 0;
  /// Derivative of the gamma function wrt inelastic strain
  virtual double dgamma(double ep, double T) const = 0;

};

/// Gamma is constant with respect to strain
class NEML_EXPORT ConstantGamma: public GammaModel {
 public:
  /// Parameter is just the constant value
  ConstantGamma(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// gamma = C
  virtual double gamma(double ep, double T) const;
  /// derivative of the gamma function
  virtual double dgamma(double ep, double T) const;

  /// Getter for the constant value
  double g(double T) const;

 private:
  const std::shared_ptr<const Interpolate> g_;
};

static Register<ConstantGamma> regConstantGamma;

/// Gamma evolves with a saturating Voce form
class NEML_EXPORT SatGamma: public GammaModel {
 public:
  /// Parameters are the initial value of gamma, the saturated value of gamma
  /// and the saturation speed constant
  SatGamma(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// gamma = gs + (g0 - gs) * exp(-beta * ep)
  virtual double gamma(double ep, double T) const;
  /// Derivative of the gamma function
  virtual double dgamma(double ep, double T) const;

  /// Parameter getter
  double gs(double T) const;
  /// Parameter getter
  double g0(double T) const;
  /// Parameter getter
  double beta(double T) const;

 private:
  const std::shared_ptr<const Interpolate> gs_, g0_, beta_;
};

static Register<SatGamma> regSatGamma;

/// Chaboche model: generalized Frederick-Armstrong
//    This model degenerates to Frederick-Armstrong for n = 1 and
//    A = 0, a = 1
class NEML_EXPORT Chaboche: public NonAssociativeHardening {
 public:
  /// Parameters: isotropic hardening model, vector of backstress constants C,
  /// vector of gamma functions, vector of static recovery constants A,
  /// vector static recovery constants a, flag trigger the nonisothermal terms
  Chaboche(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// 1 (isotropic) + 6 (backstress) = 7
  virtual size_t ninter() const; // How many "q" variables it spits out

  /// Setup the internal state
  virtual void populate_hist(History & h) const;
  /// Initialize the history
  virtual void init_hist(History & h) const;

  /// Map the isotropic variable, map the backstresses
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of map
  virtual void dq_da(const double * const alpha, double T, double * const qv) const;

  /// Hardening proportional to inelastic strain rate
  /// Assume associated isotropic hardening, each backstress is
  /// -2/3 * C * X/norm(X) - sqrt(2/3) gamma(ep) * X
  virtual void h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h wrt stress
  virtual void dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h wrt history
  virtual void dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening proportional to time
  /// Zero for the isotropic part, -A sqrt(3/2) pow(norm(X), a-1) * X
  /// for each backstress
  virtual void h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_time wrt stress
  virtual void dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_time wrt history
  virtual void dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening proportional to the temperature rate
  /// Zero for the isotropic part
  /// Zero for each backstress if noniso = False
  /// -sqrt(2/3) * dC/dT / C * X for each backstress if noniso = True
  virtual void h_temp(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_temp wrt stress
  virtual void dh_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_temp wrt history
  virtual void dh_da_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Getter for the number of backstresses
  int n() const;

  /// Getter for the C constants
  std::vector<double> c(double T) const;

 private:
  void backstress_(const double * const alpha, double * const X) const;

 private:
  std::shared_ptr<IsotropicHardeningRule> iso_;
  const std::vector<std::shared_ptr<Interpolate>> c_;  
  const size_t n_;
  std::vector<std::shared_ptr<GammaModel>> gmodels_;

  const std::vector<std::shared_ptr<Interpolate>> A_;
  const std::vector<std::shared_ptr<Interpolate>> a_;
  const bool relax_;

  const bool noniso_;
};

static Register<Chaboche> regChaboche;

/// Chaboche model with hard-coded, recoverable Voce isotropic hardening
class NEML_EXPORT ChabocheVoceRecovery: public NonAssociativeHardening {
 public:
  /// Parameters: isotropic hardening model, vector of backstress constants C,
  /// vector of gamma functions, vector of static recovery constants A,
  /// vector static recovery constants a, flag trigger the nonisothermal terms
  ChabocheVoceRecovery(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// 1 (isotropic) + 6 (backstress) = 7
  virtual size_t ninter() const; // How many "q" variables it spits out
  
  /// Setup the internal state
  virtual void populate_hist(History & h) const;
  /// Initialize the history
  virtual void init_hist(History & h) const;

  /// Map the isotropic variable, map the backstresses
  virtual void q(const double * const alpha, double T, double * const qv) const;
  /// Derivative of map
  virtual void dq_da(const double * const alpha, double T, double * const qv) const;

  /// Hardening proportional to inelastic strain rate
  /// Assume associated isotropic hardening, each backstress is
  /// -2/3 * C * X/norm(X) - sqrt(2/3) gamma(ep) * X
  virtual void h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h wrt stress
  virtual void dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h wrt history
  virtual void dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening proportional to time
  virtual void h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_time wrt stress
  virtual void dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_time wrt history
  virtual void dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening proportional to the temperature rate
  virtual void h_temp(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_temp wrt stress
  virtual void dh_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_temp wrt history
  virtual void dh_da_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Getter for the number of backstresses
  int n() const;

 private:
  void backstress_(const double * const alpha, double * const X) const;

 private:
  const std::shared_ptr<Interpolate> s0_;
  const std::shared_ptr<Interpolate> theta0_;
  const std::shared_ptr<Interpolate> Rmax_;
  const std::shared_ptr<Interpolate> Rmin_;
  const std::shared_ptr<Interpolate> r1_;
  const std::shared_ptr<Interpolate> r2_;

  const std::vector<std::shared_ptr<Interpolate>> c_;  
  const size_t n_;
  std::vector<std::shared_ptr<GammaModel>> gmodels_;

  const std::vector<std::shared_ptr<Interpolate>> A_;
  const std::vector<std::shared_ptr<Interpolate>> a_;

  const bool noniso_;
};

static Register<ChabocheVoceRecovery> regChabocheVoceRecovery;

} // namespace neml

#endif // HARDENING_H
