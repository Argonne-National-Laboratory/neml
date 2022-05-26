#ifndef RI_FLOW_H
#define RI_FLOW_H

#include "objects.h"
#include "surfaces.h"
#include "hardening.h"
#include "history.h"
#include "math/nemlmath.h"

#include "windows.h"

#include <memory>

namespace neml {

/// ABC describing rate independent flow
class NEML_EXPORT RateIndependentFlowRule: public HistoryNEMLObject {
 public:
  RateIndependentFlowRule(ParameterSet & params);
  
  /// Yield surface
  virtual void f(const double* const s, const double* const alpha, double T,
                double & fv) const = 0;
  /// Partial derivative of the surface wrt stress
  virtual void df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const = 0;
  /// Partial derivative of the surface wrt history
  virtual void df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const = 0;

  /// Flow function
  virtual void g(const double * const s, const double * const alpha, double T,
                double * const gv) const = 0;
  /// Partial derivative of the flow function wrt stress
  virtual void dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const = 0;
  /// Partial derivative of the flow function wrt history
  virtual void dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const = 0;

  /// Hardening rule
  virtual void h(const double * const s, const double * const alpha, double T,
                double * const hv) const = 0;
  /// Partial derivative of the hardening rule wrt. stress
  virtual void dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
  /// Partial derivative of the hardening rule wrt. history
  virtual void dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
};

/// Implementation of associative RI flow
class NEML_EXPORT RateIndependentAssociativeFlow: public RateIndependentFlowRule {
 public:
  /// Parameters: yield surface and hardening rule
  RateIndependentAssociativeFlow(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();

  /// Setup internal state
  virtual void populate_hist(History & h) const;
  /// Setup the history at time zero
  virtual void init_hist(History & h) const;

  /// Yield surface
  virtual void f(const double* const s, const double* const alpha, double T,
                double & fv) const;
  /// Partial derivative of the surface wrt stress
  virtual void df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const;
  /// Partial derivative of the surface wrt history
  virtual void df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const;

  /// Flow function
  virtual void g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Partial derivative of the flow function wrt stress
  virtual void dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Partial derivative of the flow function wrt history
  virtual void dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Hardening rule
  virtual void h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Partial derivative of the hardening rule wrt. stress
  virtual void dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Partial derivative of the hardening rule wrt. history
  virtual void dh_da(const double * const s, const double * const alpha, double T,
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
class NEML_EXPORT RateIndependentNonAssociativeHardening: public RateIndependentFlowRule {
 public:
  RateIndependentNonAssociativeHardening(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from a parameter set
  static ParameterSet parameters();

  /// Setup internal state
  virtual void populate_hist(History & h) const;
  /// Setup the history at time zero
  virtual void init_hist(History & h) const;

  /// Yield surface
  virtual void f(const double* const s, const double* const alpha, double T,
                double & fv) const;
  /// Partial derivative of the surface wrt stress
  virtual void df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const;
  /// Partial derivative of the surface wrt history
  virtual void df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const;

  /// Flow function
  virtual void g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  /// Partial derivative of the flow function wrt stress
  virtual void dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  /// Partial derivative of the flow function wrt history
  virtual void dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  /// Hardening rule
  virtual void h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  /// Partial derivative of the hardening rule wrt. stress
  virtual void dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  /// Partial derivative of the hardening rule wrt. history
  virtual void dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<NonAssociativeHardening> hardening_;
};

static Register<RateIndependentNonAssociativeHardening> regRateIndependentNonAssociativeHardening;

} // namespace neml

#endif // RI_FLOW_H
