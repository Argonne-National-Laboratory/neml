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
  /// Populate history with the correct variable names and types
  virtual void populate_history(History & history) const = 0;
  /// Initialize history with actual starting values
  virtual void init_history(History & history) const = 0;
  
  /// Hook to allow user to decouple parts of the update
  /// Any items added to the History object will be treated explicitly
  virtual History decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) = 0;
  
  /// Stress rate
  virtual Symmetric stress_rate(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  /// Derivative of stress rate with respect to stress
  virtual SymSymR4 d_stress_rate_d_stress(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  /// Derivative of the stress rate with respect to the deformation rate
  virtual SymSymR4 d_stress_rate_d_d(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  /// Derivative of the stress rate with respect to the vorticity
  virtual SymSkewR4 d_stress_rate_d_w(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  /// Derivative of the stress rate with respect to the history
  virtual History d_stress_rate_d_history(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  
  /// History rate
  virtual History history_rate(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  /// Derivative of the history rate with respect to the stress
  virtual History d_history_rate_d_stress(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  /// Derivative of the history rate with respect to the deformation rate
  virtual History d_history_rate_d_d(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  /// Derivative of the history rate with respect to the vorticity
  virtual History d_history_rate_d_w(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  /// Derivative of the history rate with respect to the history
  virtual History d_history_rate_d_history(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const = 0;
  
  /// Derivative of the stress rate with respect to the deformation 
  /// keeping fixed variables fixed.
  virtual SymSymR4 d_stress_rate_d_d_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed);
  /// Derivative of the stress rate with respect to the vorticity keeping
  /// fixed variables fixed
  virtual SymSkewR4 d_stress_rate_d_w_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed);
  
  /// Derivative of the history rate with respect to the deformation rate
  /// keeping fixed variables in fixed constant
  virtual History d_history_rate_d_d_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed);
  /// Derivative of the history rate with respect to the vorticity keeping
  /// fixed variables in fixed constant
  virtual History d_history_rate_d_w_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed);
  
  /// The spin rate
  virtual Skew spin(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const = 0;
  
  /// Helper to calculate elastic strains
  virtual Symmetric elastic_strains(const Symmetric & stress,
                                    const Orientation & Q,
                                    const History & history,
                                    double T) = 0;
};

/// My standard kinematic assumptions, outlined in the manual
class StandardKinematicModel: public KinematicModel {
 public:
  /// Initialize with elastic and inelastic models
  StandardKinematicModel(std::shared_ptr<LinearElasticModel> emodel,
                         std::shared_ptr<InelasticModel> imodel);
  /// Destructor
  virtual ~StandardKinematicModel();

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();
  
  /// Populate a history object with the correct variables
  virtual void populate_history(History & history) const;
  /// Initialize the history object with the starting values
  virtual void init_history(History & history) const;

  /// Anything added to this History instance will be kept fixed in the 
  /// implicit integration
  virtual History decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T);

  /// Stress rate
  virtual Symmetric stress_rate(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  /// Derivative of stress rate with respect to stress
  virtual SymSymR4 d_stress_rate_d_stress(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  /// Derivative of the stress rate with respect to the deformation rate
  virtual SymSymR4 d_stress_rate_d_d(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  /// Derivative of the stress rate with respect to the vorticity
  virtual SymSkewR4 d_stress_rate_d_w(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  /// Derivative of the stress rate with respect to the history
  virtual History d_stress_rate_d_history(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  
  /// History rate
  virtual History history_rate(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  /// Derivative of the history rate with respect to the stress
  virtual History d_history_rate_d_stress(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  /// Derivative of the history rate with respect to the deformation rate
  virtual History d_history_rate_d_d(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  /// Derivative of the history rate with respect to the vorticity
  virtual History d_history_rate_d_w(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  /// Derivative of the history rate with respect to the history
  virtual History d_history_rate_d_history(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;
  
  /// Derivative of the stress rate with respect to the vorticity keeping
  /// fixed variables fixed
  virtual SymSkewR4 d_stress_rate_d_w_decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed);
  
  /// The spin rate
  virtual Skew spin(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T) const;
  
  /// Helper to calculate elastic strains
  virtual Symmetric elastic_strains(const Symmetric & stress,
                                    const Orientation & Q,
                                    const History & history,
                                    double T);

 private:
  std::shared_ptr<LinearElasticModel> emodel_;
  std::shared_ptr<InelasticModel> imodel_;
};

static Register<StandardKinematicModel> regStandardKinematicModel;

} // namespace neml

#endif // KINEMATICS_H
