#ifndef SLIPHARDEN_H
#define SLIPHARDEN_H

#include "sliprules.h"
#include "crystallography.h"

#include "../objects.h"
#include "../history.h"
#include "../interpolate.h"

#include "../math/tensors.h"
#include "../math/rotations.h"

#include <map>
#include <vector>
#include <string>

namespace neml {

class SlipRule; // Why would we need a forward declaration?

/// ABC for a slip hardening model
class SlipHardening: public NEMLObject
{
 public:
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
                     const Lattice & L, double T, const SlipRule & R) const = 0;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress, 
                             const Orientation & Q, const History & history,
                             const Lattice & L, double T,
                             const SlipRule & R) const = 0;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress, 
                 const Orientation & Q,
                 const History & history,
                 const Lattice & L, 
                 double T, const SlipRule & R) const = 0;
};

/// Slip strength rules where all systems share the same strength
class SlipSingleHardening: public SlipHardening
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
class SlipSingleStrengthHardening: public SlipSingleHardening
{
 public:
  /// Request whatever history you will need
  virtual void populate_history(History & history) const;
  /// Setup history
  virtual void init_history(History & history) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress, 
                     const Orientation & Q, const History & history,
                     const Lattice & L, double T, const SlipRule & R) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress, 
                             const Orientation & Q, const History & history,
                             const Lattice & L, double T,
                             const SlipRule & R) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress, 
                 const Orientation & Q,
                 const History & history,
                 const Lattice & L, 
                 double T, const SlipRule & R) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T) const;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T) const;

  /// Static (not evolving) strength
  virtual double static_strength(double T) const = 0;
  
  /// Setup the scalar
  virtual double init_strength() const = 0;

  /// Scalar evolution law
  virtual double hist_rate(const Symmetric & stress, const Orientation & Q,
                           const History & history,
                           const Lattice & L, double T, const SlipRule & R) const = 0;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric & stress, const Orientation & Q, 
                                         const History & history,
                                         const Lattice & L, double T,
                                         const SlipRule & R) const = 0;
  /// Derivative of scalar law wrt the scalar
  virtual double d_hist_rate_d_strength(const Symmetric & stress, const Orientation & Q, 
                                        const History & history,
                                        const Lattice & L, double T,
                                        const SlipRule & R) const = 0;
};

/// Slip strength rule where the single strength evolves with sum|dg|
class PlasticSlipHardening: public SlipSingleStrengthHardening
{
 public:
  /// Scalar evolution law
  virtual double hist_rate(const Symmetric & stress, const Orientation & Q,
                           const History & history, const Lattice & L, double T, const SlipRule & R) const;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric & stress, const Orientation & Q, 
                                         const History & history, const Lattice & L, double T,
                                         const SlipRule & R) const;
  /// Derivative of scalar law wrt the scalar
  virtual double d_hist_rate_d_strength(const Symmetric & stress, const Orientation & Q, 
                                        const History & history, const Lattice & L, double T,
                                        const SlipRule & R) const;

  /// Prefactor
  virtual double hist_factor(double strength, const Lattice & L, double T) const = 0;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, const Lattice & L, double T) const = 0;

};

/// Everyone's favorite Voce model
class VoceSlipHardening: public PlasticSlipHardening
{
 public:
  VoceSlipHardening(std::shared_ptr<Interpolate> tau_sat,
                    std::shared_ptr<Interpolate> b,
                    std::shared_ptr<Interpolate> tau_0);

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
  virtual double hist_factor(double strength, const Lattice & L, double T) const;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, const Lattice & L, double T) const;

 private:
  std::shared_ptr<Interpolate> tau_sat_;
  std::shared_ptr<Interpolate> b_;
  std::shared_ptr<Interpolate> tau_0_;
};

static Register<VoceSlipHardening> regVoceSlipHardening;

}

#endif // SLIPHARDEN_H
