#ifndef GENERAL_FLOW
#define GENERAL_FLOW

#include "objects.h"
#include "elasticity.h"
#include "visco_flow.h"

#include "windows.h"

#include <cstddef>

namespace neml {

/// ABC for a completely general flow rule...
class NEML_EXPORT GeneralFlowRule: public HistoryNEMLObject {
 public:
  GeneralFlowRule(ParameterSet & params);

  /// Stress rate
  virtual void s(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const sdot) = 0;
  /// Partial of stress rate wrt stress
  virtual void ds_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot) = 0;
  /// Partial of stress rate wrt history
  virtual void ds_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot) = 0;
  /// Partial of stress rate wrt strain rate
  virtual void ds_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot) = 0;

  /// History rate
  virtual void a(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const adot) = 0;
  /// Partial of history rate wrt stress
  virtual void da_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot) = 0;
  /// Partial of history rate wrt history
  virtual void da_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot) = 0;
  /// Partial of history rate wrt strain rate
  virtual void da_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot) = 0;

  /// The implementation needs to define inelastic dissipation
  virtual void work_rate(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double & p_rate);

  /// The implementation needs to define elastic strain
  virtual void elastic_strains(const double * const s_np1, double T_np1,
                              double * const e_np1) const = 0;

  /// Set a new elastic model
  virtual void set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);

  /// Optional method for modifying the initial guess
  virtual void override_guess(double * const x);
};

/// Thermo-visco-plasticity
//  Stress rates, temperature rates, and time can all enter the model
class NEML_EXPORT TVPFlowRule : public GeneralFlowRule {
 public:
  /// Parameters: elastic model and a viscoplastic flow rule
  TVPFlowRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Return default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from parameter set
  static ParameterSet parameters();

  // Setup internal state
  virtual void populate_hist(History & h) const;
  /// Initialize history
  virtual void init_hist(History & h) const;

  /// Stress rate
  virtual void s(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const sdot);
  /// Partial of stress rate wrt stress
  virtual void ds_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);
  /// Partial of stress rate wrt history
  virtual void ds_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);
  /// Partial of stress rate wrt strain rate
  virtual void ds_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);

  /// History rate
  virtual void a(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const adot);
  /// Partial of history rate wrt stress
  virtual void da_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);
  /// Partial of history rate wrt history
  virtual void da_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);
  /// Partial of history rate wrt strain rate
  virtual void da_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);

  /// The implementation needs to define inelastic dissipation
  virtual void work_rate(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double & p_rate);

  /// The implementation needs to define elastic strain
  virtual void elastic_strains(const double * const s_np1, double T_np1,
                              double * const e_np1) const;

  /// Set a new elastic model
  virtual void set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);

  /// Override the initial guess
  virtual void override_guess(double * const x);

 private:
  std::shared_ptr<LinearElasticModel> elastic_;
  std::shared_ptr<ViscoPlasticFlowRule> flow_;
};

static Register<TVPFlowRule> regTVPFlowRule;

} // namespace neml

#endif
