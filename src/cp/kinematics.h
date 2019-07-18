#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "inelasticity.h"
#include "crystallography.h"

#include "../objects.h"
#include "../history.h"
#include "../elasticity.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

namespace neml {

/// Describes the stress, history, and rotation rates
class KinematicModel: public NEMLObject {
 public:
  virtual void populate_history(History & history) const = 0;
  virtual void init_history(History & history) const = 0;
  
  /// Hook to allow user to decouple parts of the update
  virtual void decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) = 0;

  virtual Symmetric stress_rate(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  virtual SymSym d_stress_rate_d_stress(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  virtual SymSym d_stress_rate_d_d(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  virtual SymSkew d_stress_rate_d_w(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  virtual History d_stress_rate_d_history(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;

  virtual History history_rate(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  virtual History d_history_rate_d_stress(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  virtual History d_history_rate_d_d(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  virtual History d_history_rate_d_w(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  virtual History d_history_rate_d_history(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;

  virtual SymSym d_stress_rate_d_d_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T);
  virtual SymSkew d_stress_rate_d_w_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T);

  virtual History d_history_rate_d_d_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T);
  virtual History d_history_rate_d_w_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T);

  virtual Skew spin(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;

  virtual Symmetric elastic_strains(const Symmetric & stress,
                                    const Orientation & Q,
                                    const History & history,
                                    double T) = 0;
};

/// My standard kinematic assumptions
class StandardKinematicModel: public KinematicModel {
 public:
  StandardKinematicModel(std::shared_ptr<LinearElasticModel> emodel,
                         std::shared_ptr<InelasticModel> imodel);
  virtual ~StandardKinematicModel();

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  virtual void populate_history(History & history) const;
  virtual void init_history(History & history) const;

  /// Hook to allow user to decouple parts of the update
  virtual void decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T);

  virtual Symmetric stress_rate(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  virtual SymSym d_stress_rate_d_stress(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  virtual SymSym d_stress_rate_d_d(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  virtual SymSkew d_stress_rate_d_w(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  virtual History d_stress_rate_d_history(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;

  virtual History history_rate(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  virtual History d_history_rate_d_stress(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  virtual History d_history_rate_d_d(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  virtual History d_history_rate_d_w(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  virtual History d_history_rate_d_history(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;

  virtual SymSkew d_stress_rate_d_w_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T);

  virtual Skew spin(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;

  virtual Symmetric elastic_strains(const Symmetric & stress,
                                    const Orientation & Q,
                                    const History & history,
                                    double T);

 private:
  std::shared_ptr<LinearElasticModel> emodel_;
  std::shared_ptr<InelasticModel> imodel_;

  Skew espin_;

  SymSym C_, S_;
};

static Register<StandardKinematicModel> regStandardKinematicModel;

} // namespace neml

#endif // KINEMATICS_H
