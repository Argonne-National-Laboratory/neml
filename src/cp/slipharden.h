#ifndef SLIPHARDEN_H
#define SLIPHARDEN_H

#include "sliprules.h"
#include "crystallography.h"

#include "../objects.h"
#include "../history.h"
#include "../interpolate.h"

#include "../math/tensors.h"
#include "../math/rotations.h"

#include "../windows.h"

#include <map>
#include <vector>
#include <string>

namespace neml {

class SlipRule; // Why would we need a forward declaration?

/// ABC for a slip hardening model
class NEML_EXPORT SlipHardening: public NEMLObject
{
 public:
  /// Report your variable names
  virtual std::vector<std::string> varnames() const = 0;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars) = 0;

  /// Request whatever history you will need
  virtual void populate_history(History & history) const = 0;
  /// Setup history
  virtual void init_history(History & history) const = 0;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             double T) const = 0;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history,
                    double T) const = 0;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R) const = 0;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R) const = 0;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R) const = 0;

 protected:
  History blank_hist() const;
};

/// Slip strength rules where all systems share the same strength
class NEML_EXPORT SlipSingleHardening: public SlipHardening
{
 public:
  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             double T) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history,
                    double T) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T) const = 0;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T) const = 0;
};

/// Slip strength rule where all systems evolve on a single scalar strength
class NEML_EXPORT SlipSingleStrengthHardening: public SlipSingleHardening
{
 public:
  SlipSingleStrengthHardening(std::string var_name = "strength");

  /// Report varnames
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// Request whatever history you will need
  virtual void populate_history(History & history) const;
  /// Setup history
  virtual void init_history(History & history) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T) const;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T) const;

  /// Set the variable's name
  void set_variable(std::string name);

  /// Static (not evolving) strength
  virtual double static_strength(double T) const = 0;

  /// Setup the scalar
  virtual double init_strength() const = 0;

  /// Scalar evolution law
  virtual double hist_rate(const Symmetric & stress, const Orientation & Q,
                           const History & history,
                           Lattice & L, double T, const SlipRule & R) const = 0;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric & stress, const Orientation & Q,
                                         const History & history,
                                         Lattice & L, double T,
                                         const SlipRule & R) const = 0;
  /// Derivative of scalar law wrt the scalar
  virtual History d_hist_rate_d_hist(const Symmetric & stress, const Orientation & Q,
                                     const History & history,
                                     Lattice & L, double T,
                                     const SlipRule & R) const = 0;

 protected:
  std::string var_name_;
};

/// Sum of individual SlipSingleStrenghHardening models (static strengths also
/// summed)
class NEML_EXPORT SumSlipSingleStrengthHardening: public SlipSingleHardening
{
 public:
  /// Initialize with a list of models
  SumSlipSingleStrengthHardening(std::vector<std::shared_ptr<SlipSingleStrengthHardening>>
                                 models);

  /// Report varnames
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Request whatever history you will need
  virtual void populate_history(History & history) const;
  /// Setup history
  virtual void init_history(History & history) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T) const;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T) const;

 private:
  size_t nmodels() const;
  const std::vector<std::shared_ptr<SlipSingleStrengthHardening>> models_;
};

static Register<SumSlipSingleStrengthHardening> regSumSlipSingleStrengthHardening;

/// Slip strength rule where the single strength evolves with sum|dg|
class NEML_EXPORT PlasticSlipHardening: public SlipSingleStrengthHardening
{
 public:
  PlasticSlipHardening(std::string var_name = "strength");

  /// Scalar evolution law
  virtual double hist_rate(const Symmetric & stress, const Orientation & Q,
                           const History & history, Lattice & L, double T, const SlipRule & R) const;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric & stress, const Orientation & Q,
                                         const History & history, Lattice & L, double T,
                                         const SlipRule & R) const;
  /// Derivative of scalar law wrt the scalar
  virtual History d_hist_rate_d_hist(const Symmetric & stress, const Orientation & Q,
                                     const History & history, Lattice & L, double T,
                                     const SlipRule & R) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice & L, double T) const = 0;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice & L, double T) const = 0;

};

/// Everyone's favorite Voce model
class NEML_EXPORT VoceSlipHardening: public PlasticSlipHardening
{
 public:
  /// Initialize with the saturated strength, the rate constant, and a constant
  /// strength
  VoceSlipHardening(std::shared_ptr<Interpolate> tau_sat,
                    std::shared_ptr<Interpolate> b,
                    std::shared_ptr<Interpolate> tau_0,
                    std::string var_name = "strength");

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup the scalar
  virtual double init_strength() const;

  /// Static strength
  virtual double static_strength(double T) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice & L, double T) const;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice & L, double T) const;

 private:
  std::shared_ptr<Interpolate> tau_sat_;
  std::shared_ptr<Interpolate> b_;
  std::shared_ptr<Interpolate> tau_0_;
};

static Register<VoceSlipHardening> regVoceSlipHardening;

}

#endif // SLIPHARDEN_H
