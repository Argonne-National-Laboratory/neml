#ifndef VISCO_FLOW_H
#define VISCO_FLOW_H

#include "objects.h"
#include "nemlerror.h"
#include "surfaces.h"
#include "hardening.h"
#include "interpolate.h"

#include "windows.h"

#include <memory>

namespace neml {

/// ABC describing viscoplastic flow
class NEML_EXPORT ViscoPlasticFlowRule: public NEMLObject {
 public:
  /// Number of history variables
  virtual size_t nhist() const = 0;
  /// Initialize history at time zero
  virtual int init_hist(double * const h) const = 0;

  /// Scalar flow rate
  virtual int y(const double* const s, const double* const alpha, double T,
                double & yv) const = 0;
  /// Derivative of scalar flow wrt stress
  virtual int dy_ds(const double* const s, const double* const alpha, double T,
                double * const dyv) const = 0;
  /// Derivative of scalar flow wrt history
  virtual int dy_da(const double* const s, const double* const alpha, double T,
                double * const dyv) const = 0;

  /// Contribution towards the flow proportional to the scalar inelastic
  /// strain rate
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const = 0;
  /// Derivative of g wrt stress
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const = 0;
  /// Derivative of g wrt history
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const = 0;

  /// Contribution towards the flow proportional directly to time
  virtual int g_time(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Derivative of g_time wrt stress
  virtual int dg_ds_time(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Derivative of g_time wrt history
  virtual int dg_da_time(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Contribution towards the flow proportional to the temperature rate
  virtual int g_temp(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Derivative of g_temp wrt stress
  virtual int dg_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Derivative of g_temp wrt history
  virtual int dg_da_temp(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Hardening rate proportional to the scalar inelastic strain rate
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const = 0;
  /// Derivative of h wrt stress
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
  /// Derivative of h wrt history
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;

  /// Hardening rate proportional directly to time
  virtual int h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_time wrt stress
  virtual int dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_time wrt history
  virtual int dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening rate proportional to the temperature rate
  virtual int h_temp(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_temp wrt. stress
  virtual int dh_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_temp wrt history
  virtual int dh_da_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
};

/// The "g" function in the Perzyna model -- often a power law
class NEML_EXPORT GFlow: public NEMLObject {
 public:
  /// The value of g
  virtual double g(double f, double T) const = 0;
  /// The derivative of g wrt to the flow surface
  virtual double dg(double f, double T) const = 0;
};

/// g is a power law
class NEML_EXPORT GPowerLaw: public GFlow {
 public:
  /// Parameter: the power law exponent
  GPowerLaw(std::shared_ptr<Interpolate> n, std::shared_ptr<Interpolate> eta);

  /// String type for the object system
  static std::string type();
  /// Default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from parameters
  static ParameterSet parameters();

  /// g = (f/eta)^n
  virtual double g(double f, double T) const;
  /// Derivative of g wrt f
  virtual double dg(double f, double T) const;

  /// Helper, just return the power law exponent
  double n(double T) const;

  /// Helper, just returns the fluidity
  double eta(double T) const;

 private:
  const std::shared_ptr<const Interpolate> n_;
  const std::shared_ptr<const Interpolate> eta_;
};

static Register<GPowerLaw> regGPowerLaw;

/// Perzyna associative viscoplasticity
class NEML_EXPORT PerzynaFlowRule : public ViscoPlasticFlowRule {
 public:
  /// Parameters: a flow surface, a hardening rule, and the rate sensitivity
  /// function
  PerzynaFlowRule(std::shared_ptr<YieldSurface> surface,
                  std::shared_ptr<HardeningRule> hardening,
                  std::shared_ptr<GFlow> g);

  /// String type for the object system
  static std::string type();
  /// Default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from parameters
  static ParameterSet parameters();

  /// Number of history variables
  virtual size_t nhist() const;
  /// Initialize history at time zero
  virtual int init_hist(double * const h) const;

  /// Scalar strain rate
  virtual int y(const double* const s, const double* const alpha, double T,
                double & yv) const;
  /// Derivative of y wrt stress
  virtual int dy_ds(const double* const s, const double* const alpha, double T,
                double * const dyv) const;
  /// Derivative of y wrt history
  virtual int dy_da(const double* const s, const double* const alpha, double T,
                double * const dyv) const;

  /// Flow rule proportional to the scalar strain rate
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Derivative of g wrt stress
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Derivative of g wrt history
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Hardening rule proportional to the scalar strain rate
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h wrt stress
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h wrt history
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<HardeningRule> hardening_;
  std::shared_ptr<GFlow> g_;
};

static Register<PerzynaFlowRule> regPerzynaFlowRule;

/// Various Chaboche type fluidity models.
//
//  These depend only on the equivalent plastic strain
class NEML_EXPORT FluidityModel: public NEMLObject {
 public:
  /// Value of viscosity as a function of temperature and inelastic strain
  virtual double eta(double a, double T) const = 0;
  /// Derivative of viscosity wrt inelastic strain
  virtual double deta(double a, double T) const = 0;
};

/// The fluidity is constant with respect to plastic strain
class NEML_EXPORT ConstantFluidity: public FluidityModel {
 public:
  /// Parameter: constant value of viscosity
  ConstantFluidity(std::shared_ptr<Interpolate> eta);

  /// String type for the object system
  static std::string type();
  /// Initialize with a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Value of eta
  virtual double eta(double a, double T) const;
  /// Derivative of eta wrt inelastic strain (zero for this implementation)
  virtual double deta(double a, double T) const;

 private:
  const std::shared_ptr<const Interpolate> eta_;
};

static Register<ConstantFluidity> regConstantFluidity;

/// Voce-like saturating fluidity
class NEML_EXPORT SaturatingFluidity: public FluidityModel {
 public:
  /// Parameters: K0, initial viscosity, A, saturated increase in viscosity,
  /// b, sets saturation rate
  SaturatingFluidity(std::shared_ptr<Interpolate> K0,
                     std::shared_ptr<Interpolate> A,
                     std::shared_ptr<Interpolate> b);

  /// String type for the object system
  static std::string type();
  /// Initialize with a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Value of eta
  virtual double eta(double a, double T) const;
  /// Derivative of eta wrt inelastic strain
  virtual double deta(double a, double T) const;

 private:
  const std::shared_ptr<const Interpolate> K0_, A_, b_;
};

static Register<SaturatingFluidity> regSaturatingFluidity;

/// Non-associative flow based on Chaboche's viscoplastic formulation
//
//  It uses an associative flow rule, a non-associative hardening rule
//  (which should be Chaboche's for the full model), and a Perzyna rate
//  rule with a potentially history-dependent fluidity
//
class NEML_EXPORT ChabocheFlowRule: public ViscoPlasticFlowRule {
 public:
  /// Parameters: a yield surface, a nonassociative hardening rule,
  /// the fluidity function, and a rate sensitivity exponent
  ChabocheFlowRule(std::shared_ptr<YieldSurface> surface,
                   std::shared_ptr<NonAssociativeHardening> hardening,
                   std::shared_ptr<FluidityModel> fluidity,
                   std::shared_ptr<Interpolate> n);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// Number of history variables (from the hardening model)
  virtual size_t nhist() const;
  /// Initialize history at time zero
  virtual int init_hist(double * const h) const;

  // Scalar inelastic strain rate
  virtual int y(const double* const s, const double* const alpha, double T,
                double & yv) const;
  /// Derivative of y wrt stress
  virtual int dy_ds(const double* const s, const double* const alpha, double T,
                double * const dyv) const;
  /// Derivative of y wrt history
  virtual int dy_da(const double* const s, const double* const alpha, double T,
                double * const dyv) const;

  /// Flow rule proportional to the scalar strain rate
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Derivative of g wrt stress
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Derivative of g wrt history
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Hardening rule proportional to the scalar strain rate
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h wrt stress
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h wrt history
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening rule proportional to time
  virtual int h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_time wrt stress
  virtual int dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_time wrt history
  virtual int dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening rule proportional to temperature rate
  virtual int h_temp(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_temp wrt stress
  virtual int dh_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_temp wrt history
  virtual int dh_da_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<NonAssociativeHardening> hardening_;
  std::shared_ptr<FluidityModel> fluidity_;
  const std::shared_ptr<const Interpolate> n_;
  const bool recovery_;
};

static Register<ChabocheFlowRule> regChabocheFlowRule;

/// Non-associative flow for Gr. 91 from Yaguchi & Takahashi (2000) + (2005)
//
//  A modified Chaboche rule with temperature interpolation even between
//  473 and 873 K.
//
//  Modifications include an evolving Chaboche constant
//
//  Interpolations are hard-coded because of their complexity
//  They are public so I can easily test them
//
class NEML_EXPORT YaguchiGr91FlowRule: public ViscoPlasticFlowRule {
 public:
  /// All parameters are hard coded to those given in the paper
  YaguchiGr91FlowRule();

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameter set
  static ParameterSet parameters();

  /// Number of history variables (14)
  virtual size_t nhist() const;
  /// Initialize history (6 values for X1, 6 values for X2, 1 value for Q and
  /// 1 value for sigma_a)
  virtual int init_hist(double * const h) const;

  /// Scalar inelastic strain rate
  virtual int y(const double* const s, const double* const alpha, double T,
                double & yv) const;
  /// Derivative of y wrt stress
  virtual int dy_ds(const double* const s, const double* const alpha, double T,
                double * const dyv) const;
  /// Derivative of y wrt history
  virtual int dy_da(const double* const s, const double* const alpha, double T,
                double * const dyv) const;

  /// Flow rule proportional to the scalar strain rate
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Derivative of g wrt stress
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Derivative of g wrt history
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Hardening rule proportional to scalar inelastic strain rate
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h wrt stress
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h wrt history
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Hardening rule proportional to time
  virtual int h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Derivative of h_time wrt stress
  virtual int dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Derivative of h_time wrt history
  virtual int dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  /// Value of parameter D
  double D(double T) const;
  /// Value of parameter n
  double n(double T) const;
  /// Value of parameter a_10
  double a10(double T) const;
  /// Value of parameter C2
  double C2(double T) const;
  /// Value of parameter a2
  double a2(double T) const;
  /// Value of parameter g1
  double g1(double T) const;
  /// Value of parameter g2
  double g2(double T) const;
  /// Value of parameter m
  double m(double T) const;
  /// Value of parameter b_r
  double br(double T) const;
  /// Value of parameter b_h
  double bh(double T) const;
  /// Value of parameter A
  double A(double T) const;
  /// Value of parameter B
  double B(double T) const;
  /// Value of parameter d
  double d(double T) const;
  /// Value of parameter q
  double q(double T) const;
  /// Value of parameter C1
  double C1(double T) const;

 private:
  // A few helpers
  double J2_(const double * const v) const;
  void dev_vec_deriv_(const double * const a, double * const b) const;

  double log_tol_ = 1.0e-15;
};

static Register<YaguchiGr91FlowRule> regYaguchiGr91FlowRule;

} // namespace neml

#endif // VISCO_FLOW_H
