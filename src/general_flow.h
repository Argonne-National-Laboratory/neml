#ifndef GENERAL_FLOW
#define GENERAL_FLOW

#include "elasticity.h"
#include "visco_flow.h"

#include <cstddef>

namespace neml {

/// ABC for a completely general flow rule...
class GeneralFlowRule {
 public:
  virtual size_t nhist() const = 0;
  virtual int init_hist(double * const h) = 0;

  virtual int s(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const sdot) = 0;
  virtual int ds_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot) = 0;
  virtual int ds_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot) = 0;
  virtual int ds_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot) = 0;

  virtual int a(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const adot) = 0;
  virtual int da_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot) = 0;
  virtual int da_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot) = 0;
  virtual int da_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot) = 0;
  
  // Needed to get plastic work
  virtual int work_rate(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double & p_rate);

  // And the elastic strains
  virtual int elastic_strains(const double * const s_np1, double T_np1,
                              double * const e_np1) const = 0; 
};

// First stab at thermo-visco-plasticity
class TVPFlowRule : public GeneralFlowRule {
 public:
  TVPFlowRule(std::shared_ptr<LinearElasticModel> elastic,
              std::shared_ptr<ViscoPlasticFlowRule> flow);

  virtual size_t nhist() const;
  virtual int init_hist(double * const h);

  virtual int s(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const sdot);
  virtual int ds_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);
  virtual int ds_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);
  virtual int ds_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);

  virtual int a(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const adot);
  virtual int da_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);
  virtual int da_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);
  virtual int da_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);

  virtual int work_rate(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double & p_dot);
  virtual int elastic_strains(const double * const s_np1, double T_np1,
                              double * const e_np1) const; 
  
 private:
  std::shared_ptr<LinearElasticModel> elastic_;
  std::shared_ptr<ViscoPlasticFlowRule> flow_;
};


} // namespace neml

#endif
