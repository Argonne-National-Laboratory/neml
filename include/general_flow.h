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
  virtual Symmetric s(const Symmetric & s, const History & alpha, 
                      const Symmetric & edot, double T, double Tdot) = 0;
  /// Partial of stress rate wrt stress
  virtual SymSymR4 ds_ds(const Symmetric & s, const History & alpha, 
                         const Symmetric & edot, double T, double Tdot) = 0;
  /// Partial of stress rate wrt history
  virtual History ds_da(const Symmetric & s, const History & alpha, 
                        const Symmetric & edot, double T, double Tdot) = 0;
  /// Partial of stress rate wrt strain rate
  virtual SymSymR4 ds_de(const Symmetric & s, const History & alpha, 
                         const Symmetric & edot, double T, double Tdot) = 0;

  /// History rate
  virtual History a(const Symmetric & s, const History & alpha, 
                    const Symmetric & edot, double T, double Tdot) = 0;
  /// Partial of history rate wrt stress
  virtual History da_ds(const Symmetric & s, const History & alpha, 
                        const Symmetric & edot, double T, double Tdot) = 0;
  /// Partial of history rate wrt history
  virtual History da_da(const Symmetric & s, const History & alpha, 
                        const Symmetric & edot, double T, double Tdot) = 0;
  /// Partial of history rate wrt strain rate
  virtual History da_de(const Symmetric & s, const History & alpha, 
                        const Symmetric & edot, double T, double Tdot) = 0;

  /// The implementation needs to define inelastic dissipation
  virtual double work_rate(const Symmetric & s, const History & alpha, 
                           const Symmetric & edot, double T, double Tdot);

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
  virtual Symmetric s(const Symmetric & s, const History & alpha, 
                      const Symmetric & edot, double T, double Tdot);
  /// Partial of stress rate wrt stress
  virtual SymSymR4 ds_ds(const Symmetric & s, const History & alpha, 
                         const Symmetric & edot, double T, double Tdot);
  /// Partial of stress rate wrt history
  virtual History ds_da(const Symmetric & s, const History & alpha, 
                        const Symmetric & edot, double T, double Tdot);
  /// Partial of stress rate wrt strain rate
  virtual SymSymR4 ds_de(const Symmetric & s, const History & alpha, 
                         const Symmetric & edot, double T, double Tdot);

  /// History rate
  virtual History a(const Symmetric & s, const History & alpha, 
                    const Symmetric & edot, double T, double Tdot);
  /// Partial of history rate wrt stress
  virtual History da_ds(const Symmetric & s, const History & alpha, 
                        const Symmetric & edot, double T, double Tdot);
  /// Partial of history rate wrt history
  virtual History da_da(const Symmetric & s, const History & alpha, 
                        const Symmetric & edot, double T, double Tdot);
  /// Partial of history rate wrt strain rate
  virtual History da_de(const Symmetric & s, const History & alpha, 
                        const Symmetric & edot, double T, double Tdot);

  /// The implementation needs to define inelastic dissipation
  virtual double work_rate(const Symmetric & s, const History & alpha, 
                           const Symmetric & edot, double T, double Tdot);

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
