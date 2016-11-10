#ifndef RI_FLOW_H
#define RI_FLOW_H

#include "surfaces.h"
#include "hardening.h"
#include "nemlmath.h"

#include <memory>

namespace neml {

/// ABC describing rate independent flow
class RateIndependentFlowRule {
 public:
  virtual size_t nhist() const = 0;
  virtual int init_hist(double * const h) const = 0;
  
  // Yield surface
  virtual int f(const double* const s, const double* const alpha, double T,
                double & fv) const = 0;
  virtual int df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const = 0;
  virtual int df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const = 0;

  // Flow rule
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const = 0;
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const = 0;
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const = 0;

  // Hardening rule
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const = 0;
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const = 0;
};

/// Implementation of associative RI flow
class RateIndependentAssociativeFlow: public RateIndependentFlowRule {
 public:
  RateIndependentAssociativeFlow(std::shared_ptr<YieldSurface> surface,
                                 std::shared_ptr<HardeningRule> hardening);

  virtual size_t nhist() const;
  virtual int init_hist(double * const h) const;
 
  virtual int f(const double* const s, const double* const alpha, double T,
                double & fv) const;
  virtual int df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const;
  virtual int df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const;

  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<HardeningRule> hardening_;

};

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
  virtual size_t nhist() const;
  virtual int init_hist(double * const h) const;
 
  virtual int f(const double* const s, const double* const alpha, double T,
                double & fv) const;
  virtual int df_ds(const double* const s, const double* const alpha, double T,
                double * const dfv) const;
  virtual int df_da(const double* const s, const double* const alpha, double T,
                double * const dfv) const;

  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<NonAssociativeHardening> hardening_;

};


} // namespace neml

#endif // RI_FLOW_H
