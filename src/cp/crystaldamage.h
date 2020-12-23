#pragma once

#include "../objects.h"
#include "../history.h"

#include "sliprules.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include "../windows.h"

namespace neml {

/// Abstract base class for slip plane damage models
class NEML_EXPORT CrystalDamageModel: public NEMLObject {
 public:
  CrystalDamageModel(std::vector<std::string> vars);

  /// Report the number of internal variables
  virtual size_t nvars() const;

  /// Report the names of the internal variables
  virtual std::vector<std::string> varnames() const;

  /// Set the internal variables to new names
  virtual void set_varnames(std::vector<std::string> names);

  /// Handle any setup that involves the lattice
  virtual void setup(Lattice & lattice);

  /// Setup whatever history variables the model requires
  virtual void populate_history(History & history) const;
  /// Initialize history
  virtual void init_history(History & history) const = 0;

  /// Returns the current projection operator
  virtual SymSymR4 projection(
      const Symmetric & stress, const History & damage, 
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T) = 0;
  /// Return the derivative of the projection operator wrt to the stress
  virtual SymSymSymR6 d_projection_d_stress(
      const Symmetric & stress, const History & damage,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T) = 0;
  /// Return the derivative of the projection operator wrt to the damage vars
  virtual History d_projection_d_history(
      const Symmetric & stress, const History & damage,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T) = 0;

  /// Damage variable rate
  virtual History damage_rate(
      const Symmetric & stress, const History & history, 
      const Orientation & Q, Lattice & lattice, 
      const SlipRule & slip, double T) const = 0;
  /// Derivative of each damage with respect to stress
  virtual History d_damage_d_stress(
      const Symmetric & stress, const History & history, 
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T) const = 0;
  /// Derivative of damage with respect to history
  virtual History d_damage_d_history(
      const Symmetric & stress, const History & history,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T) const = 0;

 protected:
  std::vector<std::string> varnames_;
};

/// Temp class to check interface (delete later)
class NilDamageModel: public CrystalDamageModel {
 public:
  NilDamageModel();

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Initialize history
  virtual void init_history(History & history) const;

  /// Returns the current projection operator
  virtual SymSymR4 projection(
      const Symmetric & stress, const History & damage, 
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T);
  /// Return the derivative of the projection operator wrt to the stress
  virtual SymSymSymR6 d_projection_d_stress(
      const Symmetric & stress, const History & damage,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T);
  /// Return the derivative of the projection operator wrt to the damage vars
  virtual History d_projection_d_history(
      const Symmetric & stress, const History & damage,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T);

  /// Damage along each slip plane
  virtual History damage_rate(
      const Symmetric & stress, const History & history, 
      const Orientation & Q, Lattice & lattice, 
      const SlipRule & slip, double T) const;
  /// Derivative of each damage with respect to stress
  virtual History d_damage_d_stress(
      const Symmetric & stress, const History & history, 
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T) const;
  /// Derivative of damage with respect to history
  virtual History d_damage_d_history(
      const Symmetric & stress, const History & history,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T) const;
};

static Register<NilDamageModel> regNilDamageModel;

} // namespace neml
