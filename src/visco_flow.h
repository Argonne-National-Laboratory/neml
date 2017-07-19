#ifndef VISCO_FLOW_H
#define VISCO_FLOW_H

#include "nemlerror.h"
#include "surfaces.h"
#include "hardening.h"
#include "interpolate.h"

#include <memory>

namespace neml {

/// The "g" function in the Perzyna model -- often a power law
class GFlow {
 public:
  virtual double g(double f, double T) const = 0;
  virtual double dg(double f, double T) const = 0;
};

/// g is a power law
class GPowerLaw: public GFlow {
 public:
  GPowerLaw(double n);
  GPowerLaw(std::shared_ptr<Interpolate> n);
  virtual double g(double f, double T) const;
  virtual double dg(double f, double T) const;

  double n(double T) const;

 private:
  const std::shared_ptr<const Interpolate> n_;

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

  // Flow rule wrt to strain rate
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const = 0;
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const = 0;
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const = 0;
  
  // Flow rule wrt to time
  virtual int g_time(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  virtual int dg_ds_time(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  virtual int dg_da_time(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  // Flow rule wrt to temperature
  virtual int g_temp(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  virtual int dg_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  virtual int dg_da_temp(const double * const s, const double * const alpha, double T,
               double * const dgv) const;

  // Hardening rule wrt to strain rate
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

/// Perzyna associative viscoplasticity
class PerzynaFlowRule : public ViscoPlasticFlowRule {
 public:
  PerzynaFlowRule(std::shared_ptr<YieldSurface> surface,
                  std::shared_ptr<HardeningRule> hardening,
                  std::shared_ptr<GFlow> g,
                  double eta);
  PerzynaFlowRule(std::shared_ptr<YieldSurface> surface,
                  std::shared_ptr<HardeningRule> hardening,
                  std::shared_ptr<GFlow> g,
                  std::shared_ptr<Interpolate> eta);

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
  double eta(double T) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<HardeningRule> hardening_;
  std::shared_ptr<GFlow> g_;
  const std::shared_ptr<const Interpolate> eta_;

};

/// Various Chaboche type fluidity models.
//
//  These depend only on the equivalent plastic strain
class FluidityModel {
 public:
  virtual double eta(double a, double T) const = 0;
  virtual double deta(double a, double T) const = 0;
};

class ConstantFluidity: public FluidityModel {
 public:
  ConstantFluidity(double eta);
  ConstantFluidity(std::shared_ptr<Interpolate> eta);
  virtual double eta(double a, double T) const;
  virtual double deta(double a, double T) const;

 private:
  const std::shared_ptr<const Interpolate> eta_;

};

class SaturatingFluidity: public FluidityModel {
 public:
  SaturatingFluidity(double K0, double A, double b);
  SaturatingFluidity(std::shared_ptr<Interpolate> K0,
                     std::shared_ptr<Interpolate> A,
                     std::shared_ptr<Interpolate> b);
  virtual double eta(double a, double T) const;
  virtual double deta(double a, double T) const;

 private:
  const std::shared_ptr<const Interpolate> K0_, A_, b_;

};

/// Non-associative flow based on Chaboche's viscoplastic formulation
//
//  It uses an associative flow rule, a non-associative hardening rule 
//  (which should be Chaboche's for the full model), and a Perzyna rate
//  rule with a potentially history-dependent fluidity
//
class ChabocheFlowRule: public ViscoPlasticFlowRule {
 public:
  // No recovery
  ChabocheFlowRule(std::shared_ptr<YieldSurface> surface,
                   std::shared_ptr<NonAssociativeHardening> hardening,
                   std::shared_ptr<FluidityModel> fluidity,
                   double n);
  ChabocheFlowRule(std::shared_ptr<YieldSurface> surface,
                   std::shared_ptr<NonAssociativeHardening> hardening,
                   std::shared_ptr<FluidityModel> fluidity,
                   std::shared_ptr<Interpolate> n);

  // Recovery

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

  // Hardening rule wrt to time
  virtual int h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<NonAssociativeHardening> hardening_;
  std::shared_ptr<FluidityModel> fluidity_;
  const std::shared_ptr<const Interpolate> n_;
  const bool recovery_;

};

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
class YaguchiGr91FlowRule: public ViscoPlasticFlowRule {
 public:
  YaguchiGr91FlowRule();

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

  // Hardening rule wrt to time
  virtual int h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual int dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual int dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  
  // Huge number of temperature-dependent constants
  double D(double T) const;
  double n(double T) const;
  double a10(double T) const;
  double C2(double T) const;
  double a2(double T) const;
  double g1(double T) const;
  double g2(double T) const;
  double m(double T) const;
  double br(double T) const;
  double bh(double T) const;
  double A(double T) const;
  double B(double T) const;
  double d(double T) const;
  double q(double T) const;
  double C1(double T) const;

 private:
  // A few helpers
  double J2_(const double * const v) const;
  void dev_vec_deriv_(const double * const a, double * const b) const;

  double log_tol_ = 1.0e-15;

};

} // namespace neml

#endif // VISCO_FLOW_H
