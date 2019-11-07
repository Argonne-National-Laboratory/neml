#ifndef RI_FLOW_H
#define RI_FLOW_H

#include "objects.h"
#include "surfaces.h"
#include "hardening.h"
#include "math/nemlmath.h"

#include "windows.h"

#include <memory>

namespace neml {

/// ABC describing rate independent flow
class RateIndependentFlowRule: public NEMLObject {
 public:
  /// Number of history variables
  virtual size_t nhist() const = 0;
  /// Setup the history at time zero
  virtual int init_hist(double * const h) const = 0;

  /// Yield surface
  virtual int f(const double* const s, const double* const alpha, double T,
                double & fv) const = 0;
  /// Partial derivative of the surface wrt stress
  virtual int df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const = 0;
  /// Partial derivative of the surface wrt history
  virtual int df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const = 0;

  /// Flow function
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const = 0;
  /// Partial derivative of the flow function wrt stress
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const = 0;
  /// Partial derivative of the flow function wrt history
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const = 0;

  /// Hardening rule
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const = 0;
  /// Partial derivative of the hardening rule wrt. stress
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
  /// Partial derivative of the hardening rule wrt. history
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
};

/// Implementation of associative RI flow
class RateIndependentAssociativeFlow: public RateIndependentFlowRule {
 public:
  /// Parameters: yield surface and hardening rule
  RateIndependentAssociativeFlow(std::shared_ptr<YieldSurface> surface,
                                 std::shared_ptr<HardeningRule> hardening);

  /// String type for the object system
  static std::string type();
  /// Setup default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();

  /// History size according to the HardeningRule
  virtual size_t nhist() const;
  /// Initialize history with the HardeningRule
  virtual int init_hist(double * const h) const;

  /// Yield surface
  virtual int f(const double* const s, const double* const alpha, double T,
                double & fv) const;
  /// Partial derivative of the surface wrt stress
  virtual int df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const;
  /// Partial derivative of the surface wrt history
  virtual int df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const;

  /// Flow function
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Partial derivative of the flow function wrt stress
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Partial derivative of the flow function wrt history
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Hardening rule
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Partial derivative of the hardening rule wrt. stress
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Partial derivative of the hardening rule wrt. history
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<HardeningRule> hardening_;
};

static Register<RateIndependentAssociativeFlow> regRateIndependentAssociativeFlow;

/// Associative plastic flow but non-associative hardening
//    This is a general and fairly common case where the plastic flow rule
//    (plastic strain rate evolution equation) is associative (normal) with
//    some yield surface but the hardening rule is general.
//
//    Examples of this kind of model include Frederick-Armstrong and the
//    Chaboche rate-independent varieties.
//
class RateIndependentNonAssociativeHardening: public RateIndependentFlowRule {
 public:
  RateIndependentNonAssociativeHardening(std::shared_ptr<YieldSurface> surface,
                                         std::shared_ptr<NonAssociativeHardening> hardening);

  /// String type for the object system
  static std::string type();
  /// Setup default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();

  /// History size according to the HardeningRule
  virtual size_t nhist() const;
  /// Initialize history with the HardeningRule
  virtual int init_hist(double * const h) const;

  /// Yield surface
  virtual int f(const double* const s, const double* const alpha, double T,
                double & fv) const;
  /// Partial derivative of the surface wrt stress
  virtual int df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const;
  /// Partial derivative of the surface wrt history
  virtual int df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const;

  /// Flow function
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Partial derivative of the flow function wrt stress
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Partial derivative of the flow function wrt history
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Hardening rule
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Partial derivative of the hardening rule wrt. stress
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Partial derivative of the hardening rule wrt. history
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<NonAssociativeHardening> hardening_;
};

static Register<RateIndependentNonAssociativeHardening> regRateIndependentNonAssociativeHardening;

} // namespace neml

#endif // RI_FLOW_H
