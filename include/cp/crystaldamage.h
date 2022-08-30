#pragma once

#include "../objects.h"
#include "../history.h"

#include "sliprules.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include "../windows.h"

namespace neml {

// Forward declarations
class SlipPlaneDamage;
class TransformationFunction;

/// Abstract base class for slip plane damage models
class NEML_EXPORT CrystalDamageModel: public HistoryNEMLObject {
 public:
  CrystalDamageModel(ParameterSet & params,
                     std::vector<std::string> vars);

  /// Report the number of internal variables
  virtual size_t nvars() const;

  /// Report the names of the internal variables
  virtual std::vector<std::string> varnames() const;

  /// Set the internal variables to new names
  virtual void set_varnames(std::vector<std::string> names);

  /// Setup whatever history variables the model requires
  virtual void populate_hist(History & history) const;
  /// Initialize history
  virtual void init_hist(History & history) const = 0;

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
      const SlipRule & slip, double T, const History & fixed) const = 0;
  /// Derivative of each damage with respect to stress
  virtual History d_damage_d_stress(
      const Symmetric & stress, const History & history, 
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T, const History & fixed) const = 0;
  /// Derivative of damage with respect to history
  virtual History d_damage_d_history(
      const Symmetric & stress, const History & history,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T, const History & fixed) const = 0;

 protected:
  std::vector<std::string> varnames_;
};

/// Temp class to check interface (delete later)
class NEML_EXPORT NilDamageModel: public CrystalDamageModel {
 public:
  NilDamageModel(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Initialize history
  virtual void init_hist(History & history) const;

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
      const SlipRule & slip, double T, const History & fixed) const;
  /// Derivative of each damage with respect to stress
  virtual History d_damage_d_stress(
      const Symmetric & stress, const History & history, 
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T, const History & fixed) const;
  /// Derivative of damage with respect to history
  virtual History d_damage_d_history(
      const Symmetric & stress, const History & history,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T, const History & fixed) const;
};

static Register<NilDamageModel> regNilDamageModel;

/// Project damage on each plane proportional to some damage measure on the
/// plane
class NEML_EXPORT PlanarDamageModel: public CrystalDamageModel {
 public:
  PlanarDamageModel(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Initialize history
  virtual void init_hist(History & history) const;

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
      const SlipRule & slip, double T, const History & fixed) const;
  /// Derivative of each damage with respect to stress
  virtual History d_damage_d_stress(
      const Symmetric & stress, const History & history, 
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T, const History & fixed) const;
  /// Derivative of damage with respect to history
  virtual History d_damage_d_history(
      const Symmetric & stress, const History & history,
      const Orientation & Q, Lattice & lattice,
      const SlipRule & slip, double T, const History & fixed) const;

 protected:
  History damage_history_(const History & total) const;
  History inelastic_history_(const History & total) const;


 private:
  std::shared_ptr<SlipPlaneDamage> damage_;
  std::shared_ptr<TransformationFunction> shear_transform_;
  std::shared_ptr<TransformationFunction> normal_transform_;
  std::shared_ptr<Lattice> lattice_;
};

static Register<PlanarDamageModel> regPlanarDamageModel;

/// Slip plane damage functions
class NEML_EXPORT SlipPlaneDamage : public NEMLObject {
 public:
  SlipPlaneDamage(ParameterSet & params);

  /// Initial value
  virtual double setup() const = 0;

  /// Damage rate
  virtual double damage_rate(const std::vector<double> & shears,
                             const std::vector<double> & sliprates,
                             double normal_stress,
                             double damage) = 0;
  /// Derivative wrt shears
  virtual std::vector<double> d_damage_rate_d_shear(
      const std::vector<double> & shears, const std::vector<double> & sliprates,
      double normal_stress, double damage) = 0;

  /// Derivative wrt slip rates
  virtual std::vector<double> d_damage_rate_d_slip(
      const std::vector<double> & shears, const std::vector<double> & sliprates,
      double normal_stress, double damage) = 0;

  /// Derivative wrt the normal stress
  virtual double d_damage_rate_d_normal(
      const std::vector<double> & shears, const std::vector<double> & sliprates,
      double normal_stress, double damage) = 0;

  /// Derivative wrt the damage variable
  virtual double d_damage_rate_d_damage(
      const std::vector<double> & shears, const std::vector<double> & sliprates,
      double normal_stress, double damage) = 0;
};

/// Accumulated work
class NEML_EXPORT WorkPlaneDamage : public SlipPlaneDamage
{
 public:
  WorkPlaneDamage(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Initial value
  virtual double setup() const;

  /// Damage rate
  virtual double damage_rate(const std::vector<double> & shears,
                             const std::vector<double> & sliprates,
                             double normal_stress,
                             double damage);
  /// Derivative wrt shears
  virtual std::vector<double> d_damage_rate_d_shear(
      const std::vector<double> & shears, const std::vector<double> & sliprates,
      double normal_stress, double damage);

  /// Derivative wrt slip rates
  virtual std::vector<double> d_damage_rate_d_slip(
      const std::vector<double> & shears, const std::vector<double> & sliprates,
      double normal_stress, double damage);

  /// Derivative wrt the normal stress
  virtual double d_damage_rate_d_normal(
      const std::vector<double> & shears, const std::vector<double> & sliprates,
      double normal_stress, double damage);

  /// Derivative wrt the damage variable
  virtual double d_damage_rate_d_damage(
      const std::vector<double> & shears, const std::vector<double> & sliprates,
      double normal_stress, double damage);
};

static Register<WorkPlaneDamage> regWorkPlaneDamage;

/// Transformation functions: map the damage variable + ancillary info into 
/// the range [0,1]
class NEML_EXPORT TransformationFunction: public NEMLObject {
 public:
  TransformationFunction(ParameterSet & params);

  /// Map from damage and the normal stress to [0,1]
  virtual double map(double damage, double normal_stress) = 0;
  /// Derivative of the map with respect to the damage
  virtual double d_map_d_damage(double damage, double normal_stress) = 0;
  /// Derivative of the map with respect to the normal stress
  virtual double d_map_d_normal(double damage, double normal_stress) = 0;
};

/// Sigmoid function.  x=0 -> y=0, x=c -> y=1, beta controls smoothing
class NEML_EXPORT SigmoidTransformation: public TransformationFunction {
 public:
  SigmoidTransformation(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();
  
  /// Map from damage and the normal stress to [0,1]
  virtual double map(double damage, double normal_stress);
  /// Derivative of the map with respect to damage
  virtual double d_map_d_damage(double damage, double normal_stress);
  /// Derivative of the map with respect to the normal stress
  virtual double d_map_d_normal(double damage, double normal_stress);

 private:
  double c_;
  double beta_;
  double cut_;
};

static Register<SigmoidTransformation> regSigmoidTransformation;

/// Normal stress switch: don't damage compression
class NEML_EXPORT SwitchTransformation: public TransformationFunction {
 public:
  SwitchTransformation(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Map from damage and the normal stress to [0,1]  
  virtual double map(double damage, double normal_stress);
  /// Derivative of the map with respect to damage
  virtual double d_map_d_damage(double damage, double normal_stress);
  /// Derivative of the map with respect to the normal stress
  virtual double d_map_d_normal(double damage, double normal_stress);

 private:
  std::shared_ptr<TransformationFunction> base_;
};

static Register<SwitchTransformation> regSwitchTransformation;

} // namespace neml
