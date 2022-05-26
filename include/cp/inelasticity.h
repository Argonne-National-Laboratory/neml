#ifndef INELASTICITY_H
#define INELASTICITY_H

#include "crystallography.h"
#include "sliprules.h"

#include "../objects.h"
#include "../history.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include "../windows.h"

namespace neml {

/// A inelastic model supplying D_p and W_p
//  A complete model for implicit rotations would include the derivatives wrt
//  the orientation.  Similarly it would include all the w_p derivatives.
class NEML_EXPORT InelasticModel: public HistoryNEMLObject {
 public:
  InelasticModel(ParameterSet & params);

  /// Helper for external models that want an average strength
  virtual double strength(const History & history, Lattice & L, double T,
                          const History & fixed) const = 0;

  /// Symmetric part of the plastic deformation
  virtual Symmetric d_p(const Symmetric & stress, const Orientation & Q,
                        const History & history,
                        Lattice & lattice, double T,
                        const History & fixed) const = 0;
  /// Derivative of the symmetric part with respect to stress
  virtual SymSymR4 d_d_p_d_stress(const Symmetric & stress, const Orientation & Q,
                                const History & history,
                                Lattice & lattice, double T,
                                const History & fixed) const = 0;
  /// Derivative of the symmetric part with respect to stress
  virtual History d_d_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice, double T,
                                  const History & fixed) const = 0;

  /// History rate
  virtual History history_rate(const Symmetric & stress, const Orientation & Q,
                               const History & history,
                               Lattice & lattice, double T,
                               const History & fixed) const = 0;
  /// Derivative of the history rate with respect to stress
  virtual History d_history_rate_d_stress(const Symmetric & stress,
                                          const Orientation & Q,
                                          const History & history,
                                          Lattice & lattice,
                                          double T, const History & fixed) const = 0;
  /// Derivative of the history rate with respect to the history
  virtual History d_history_rate_d_history(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & lattice,
                                         double T, const History & fixed) const = 0;
  /// Skew part of the plastic deformation rate
  virtual Skew w_p(const Symmetric & stress,
                   const Orientation & Q,
                   const History & history,
                   Lattice & lattice,
                   double T, const History & fixed) const = 0;
  /// Derivative of the skew part with respect to stress
  virtual SkewSymR4 d_w_p_d_stress(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & lattice,
                                 double T, const History & fixed) const = 0;
  /// Derivative of the skew part with respect to history
  virtual History d_w_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T, const History & fixed) const = 0;

  /// Whether this model uses the nye tensor
  virtual bool use_nye() const;
};

/// This model returns zero for the plastic deformation, resulting model
/// would be linear-elastic
class NEML_EXPORT NoInelasticity: public InelasticModel {
 public:
  /// Don't need any parameters to return zero!
  NoInelasticity(ParameterSet & params);
  /// Destructor
  virtual ~NoInelasticity();

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Add history variables (none needed)
  virtual void populate_hist(History & history) const;
  /// Define initial history (none)
  virtual void init_hist(History & history) const;

  /// Helper for external models that want an average strength
  virtual double strength(const History & history, Lattice & L, double T,
                          const History & fixed) const;

  /// Symmetric part of the plastic deformation = 0
  virtual Symmetric d_p(const Symmetric & stress,
                        const Orientation & Q,
                        const History & history,
                        Lattice & lattice,
                        double T, const History & fixed) const;
  /// Derivative of the symmetric part with respect to stress (=0)
  virtual SymSymR4 d_d_p_d_stress(const Symmetric & stress,
                                const Orientation & Q,
                                const History & history,
                                Lattice & lattice,
                                double T, const History & fixed) const;
  /// Derivative of the symmetric part with respect to history (null)
  virtual History d_d_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T, const History & fixed) const;

  /// History rate (null, as there are no history variables)
  virtual History history_rate(const Symmetric & stress, const Orientation & Q,
                               const History & history,
                               Lattice & lattice, double T,
                               const History & fixed) const;
  /// Derivative of the history rate with respect to stress (null)
  virtual History d_history_rate_d_stress(const Symmetric & stress,
                                          const Orientation & Q,
                                          const History & history,
                                          Lattice & lattice,
                                          double T, const History & fixed) const;
  /// Derivative of the history rate with respect to history (null)
  virtual History d_history_rate_d_history(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & lattice,
                                         double T, const History & fixed) const;
  // Skew part of the plastic deformation (=0)
  virtual Skew w_p(const Symmetric & stress,
                   const Orientation & Q,
                   const History & history,
                   Lattice & lattice, double T,
                   const History & fixed) const;
  /// Derivative of the skew part with respect to stress (=0)
  virtual SkewSymR4 d_w_p_d_stress(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & lattice,
                                 double T, const History & fixed) const;
  /// Derivative of the skew part with respect to history (null)
  virtual History d_w_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T, const History & fixed) const;
};

static Register<NoInelasticity> regNoInelasticity;

/// The classic crystal plasticity inelastic model, as described in the manual
class NEML_EXPORT AsaroInelasticity: public InelasticModel {
 public:
  /// Provide a SlipRule defining the slip rate/strength relations
  AsaroInelasticity(ParameterSet & params);
  /// Destructor
  virtual ~AsaroInelasticity();

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Populate the history, deferred to the SlipRule
  virtual void populate_hist(History & history) const;
  /// Initialize the history with the starting values, deferred to the SlipRule
  virtual void init_hist(History & history) const;

  /// Helper for external models that want an average strength
  virtual double strength(const History & history, Lattice & L, double T,
                          const History & fixed) const;

  /// Symmetric part of the plastic deformation rate
  virtual Symmetric d_p(const Symmetric & stress,
                        const Orientation & Q,
                        const History & history,
                        Lattice & lattice, double T,
                        const History & fixed) const;
  /// Derivative of the symmetric part with respect to stress
  virtual SymSymR4 d_d_p_d_stress(const Symmetric & stress,
                                const Orientation & Q,
                                const History & history,
                                Lattice & lattice, double T,
                                const History & fixed) const;
  /// Derivative of the symmetric part with respect to history
  virtual History d_d_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice, double T, 
                                  const History & fixed) const;

  /// History rate, deferred to the SlipRule
  virtual History history_rate(const Symmetric & stress, const Orientation & Q,
                               const History & history,
                               Lattice & lattice, double T,
                               const History & fixed) const;
  /// Derivative of the history rate with respect to stress
  virtual History d_history_rate_d_stress(const Symmetric & stress,
                                          const Orientation & Q,
                                          const History & history,
                                          Lattice & lattice, double T,
                                          const History & fixed) const;
  /// Derivative of the history rate with respect to history
  virtual History d_history_rate_d_history(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & lattice, double T,
                                         const History & fixed) const;

  /// Skew part of the plastic deformation rate
  virtual Skew w_p(const Symmetric & stress,
                   const Orientation & Q,
                   const History & history,
                   Lattice & lattice, double T,
                   const History & fixed) const;
  /// Derivative of the skew part with respect to stress
  virtual SkewSymR4 d_w_p_d_stress(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & lattice,
                                 double T,
                                 const History & fixed) const;
  /// Derivative of the skew part with respect to history
  virtual History d_w_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T,
                                  const History & fixed) const;

  /// Whether this model uses the Nye tensor
  virtual bool use_nye() const;

  /// Access to the slip rule for other models to get detailed slip information
  const SlipRule & slip_rule() const {return *rule_;};

 private:
  std::shared_ptr<SlipRule> rule_;
};

static Register<AsaroInelasticity> regAsaroInelasticity;

/// An isotropic power law creep model,
/// typically combined with slip system models to represent diffusion
class NEML_EXPORT PowerLawInelasticity: public InelasticModel {
 public:
  PowerLawInelasticity(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup history variables (none used in this implementation)
  virtual void populate_hist(History & history) const;
  /// Initialize the history variables (n/a)
  virtual void init_hist(History & history) const;

  /// Helper for external models that want an average strength
  virtual double strength(const History & history, Lattice & L, double T,
                          const History & fixed) const;

  /// Symmetric part of the deformation rate
  virtual Symmetric d_p(const Symmetric & stress,
                        const Orientation & Q,
                        const History & history,
                        Lattice & lattice,
                        double T,
                        const History & fixed) const;
  /// Derivative of the symmetric part with respect to the stress
  virtual SymSymR4 d_d_p_d_stress(const Symmetric & stress,
                                const Orientation & Q,
                                const History & history,
                                Lattice & lattice,
                                double T, const History & fixed) const;
  /// Derivative of the symmetric part with respect to the history (null here)
  virtual History d_d_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T,
                                  const History & fixed) const;

  /// History rate (null)
  virtual History history_rate(const Symmetric & stress, const Orientation & Q,
                               const History & history,
                               Lattice & lattice, double T,
                               const History & fixed) const;
  /// Derivative of the history with respect to the stress (null here)
  virtual History d_history_rate_d_stress(const Symmetric & stress,
                                          const Orientation & Q,
                                          const History & history,
                                          Lattice & lattice,
                                          double T,
                                          const History & fixed) const;
  /// Derivative of the history rate with respect to the history (null here)
  virtual History d_history_rate_d_history(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & lattice,
                                         double T,
                                         const History & fixed) const;

  /// Skew part of the deformation rate (zero here)
  virtual Skew w_p(const Symmetric & stress,
                   const Orientation & Q,
                   const History & history,
                   Lattice & lattice, double T,
                   const History & fixed) const;
  /// Derivative of the skew part with respect to stress (=0)
  virtual SkewSymR4 d_w_p_d_stress(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & lattice,
                                 double T,
                                 const History & fixed) const;
  /// Derivative of the skew part with respect to the history (null here)
  virtual History d_w_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T,
                                  const History & fixed) const;

 private:
  double seq_(const Symmetric & stress) const;

 private:
  std::shared_ptr<Interpolate> A_;
  std::shared_ptr<Interpolate> n_;
};

static Register<PowerLawInelasticity> regPowerLawInelasticity;

/// Metamodel that combines the rates of several individual InelasticModels
// The model plastic deformation rates are summed, the model histories are
// concatenated
class NEML_EXPORT CombinedInelasticity: public InelasticModel {
 public:
  /// Initialize with the list of models
  CombinedInelasticity(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Helper for external models that want an average strength
  virtual double strength(const History & history, Lattice & L, double T,
                          const History & fixed) const;

  /// Setup all history variables
  virtual void populate_hist(History & history) const;
  /// Initialize history with actual values
  virtual void init_hist(History & history) const;

  /// Sum the symmetric parts of the plastic deformation rates
  virtual Symmetric d_p(const Symmetric & stress,
                        const Orientation & Q,
                        const History & history,
                        Lattice & lattice,
                        double T, const History & fixed) const;
  /// Sum the derivatives of the symmetric part with respect to stress
  virtual SymSymR4 d_d_p_d_stress(const Symmetric & stress,
                                const Orientation & Q,
                                const History & history,
                                Lattice & lattice,
                                double T, const History & fixed) const;
  /// Concatenate the derivatives of the symmetric part with respect to history
  virtual History d_d_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T, const History & fixed) const;

  /// Concatenate the history rates
  virtual History history_rate(const Symmetric & stress, const Orientation & Q,
                               const History & history,
                               Lattice & lattice, double T, 
                               const History & fixed) const;
  /// Concatenate the derivative of the history rates with respect to stress
  virtual History d_history_rate_d_stress(const Symmetric & stress,
                                          const Orientation & Q,
                                          const History & history,
                                          Lattice & lattice,
                                          double T,
                                          const History & fixed) const;
  /// Concatenate the derivative of the history rates with respect to history
  virtual History d_history_rate_d_history(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & lattice,
                                         double T,
                                         const History & fixed) const;

  /// Sum the skew parts of the plastic deformation rates
  virtual Skew w_p(const Symmetric & stress,
                   const Orientation & Q,
                   const History & history,
                   Lattice & lattice, double T,
                   const History & fixed) const;
  /// Sum the derivatives of the skew parts with respect to the stress
  virtual SkewSymR4 d_w_p_d_stress(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & lattice,
                                 double T, const History & fixed) const;
  /// Concatenate the derivatives of the skew parts with respect to the history
  virtual History d_w_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T, const History & fixed) const;

  /// Whether this model uses the Nye tensor
  virtual bool use_nye() const;

 private:
  std::vector<std::shared_ptr<InelasticModel>> models_;
};

static Register<CombinedInelasticity> regCombinedInelasticity;

} // namespace neml

#endif // INELASTICITY_H
