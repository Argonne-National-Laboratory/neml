#ifndef VISCO_FLOW_H
#define VISCO_FLOW_H

#include "nemlerror.h"
#include "surfaces.h"
#include "hardening.h"

#include <memory>

namespace neml {

/// The "g" function in the Perzyna model -- often a power law
class GFlow {
 public:
  virtual double g(double f) const = 0;
  virtual double dg(double f) const = 0;
};

class GPowerLaw {
 public:
  GPowerLaw(double n);
  virtual double g(double f) const;
  virtual double dg(double f) const;

  double n() const;

 private:
  const double n_;

};


/// ABC describing viscoplastic flow
class ViscoPlasticFlowRule {
 public:
  virtual size_t nhist() const = 0;
  virtual int init_hist(double * const h) const = 0;
  
  // Rate rule
  virtual int y(const double* const s, const double* const alpha, double T,
                double & yv) const = 0;
  virtual int dy_ds(const double* const s, const double* const alpha, double T,
                double * const dyv) const = 0;
  virtual int dy_da(const double* const s, const double* const alpha, double T,
                double * const dyv) const = 0;

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

/// Perzyna associative viscoplasticity
class PerzynaFlowRule : public ViscoPlasticFlowRule {
 public:
  PerzynaFlowRule(std::shared_ptr<YieldSurface> surface,
                  std::shared_ptr<HardeningRule> hardening,
                  std::shared_ptr<GFlow> g,
                  double eta);

  virtual size_t nhist() const;
  virtual int init_hist(double * const h) const;
  
  // Rate rule
  virtual int y(const double* const s, const double* const alpha, double T,
                double & yv) const;
  virtual int dy_ds(const double* const s, const double* const alpha, double T,
                double * const dyv) const;
  virtual int dy_da(const double* const s, const double* const alpha, double T,
                double * const dyv) const;

  // Flow rule
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  // Hardening rule
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

  // Getter
  double eta() const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<HardeningRule> hardening_;
  std::shared_ptr<GFlow> g_;
  const double eta_;

};

/// Various Chaboche type fluidity models.
//
//  These depend only on the equivalent plastic strain
class FluidityModel {
 public:
  virtual double eta(double a) const = 0;
  virtual double deta(double a) const = 0;
};

class ConstantFluidity: public FluidityModel {
 public:
  ConstantFluidity(double eta);
  virtual double eta(double a) const;
  virtual double deta(double a) const;

 private:
  const double eta_;

};

/// Non-associative flow based on Chaboche's viscoplastic formulation
//
//  It uses an associative flow rule, a non-associative hardening rule 
//  (which should be Chaboche's for the full model), and a Perzyna rate
//  rule with a potentially history-dependent fluidity
//
class ChabocheFlowRule: public ViscoPlasticFlowRule {
 public:
  ChabocheFlowRule(std::shared_ptr<YieldSurface> surface,
                   std::shared_ptr<NonAssociativeHardening> hardening,
                   std::shared_ptr<FluidityModel> fluidity,
                   double n);

  virtual size_t nhist() const;
  virtual int init_hist(double * const h) const;
  
  // Rate rule
  virtual int y(const double* const s, const double* const alpha, double T,
                double & yv) const;
  virtual int dy_ds(const double* const s, const double* const alpha, double T,
                double * const dyv) const;
  virtual int dy_da(const double* const s, const double* const alpha, double T,
                double * const dyv) const;

  // Flow rule
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  // Hardening rule
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<NonAssociativeHardening> hardening_;
  std::shared_ptr<FluidityModel> fluidity_;
  const double n_;

};

} // namespace neml

#endif // VISCO_FLOW_H
