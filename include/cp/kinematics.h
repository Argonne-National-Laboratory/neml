#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "inelasticity.h"
#include "crystallography.h"
#include "crystaldamage.h"

#include "../objects.h"
#include "../history.h"
#include "../elasticity.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include "../windows.h"

namespace neml {

/// Describes the stress, history, and rotation rates
class NEML_EXPORT KinematicModel: public HistoryNEMLObject {
 public:
  KinematicModel(ParameterSet & params);

  /// Helper for external models that want a strength
  virtual double strength(const History & history, Lattice & L, double T,
                          const History & fixed) const = 0;

  /// Hook to allow user to decouple parts of the update
  /// Any items added to the History object will be treated explicitly
  virtual History decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) = 0;

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
      double T, const History & fixed) const = 0;

  /// Helper to calculate elastic strains
  virtual Symmetric elastic_strains(const Symmetric & stress,
                                    Lattice & lattice,
                                    const Orientation & Q,
                                    const History & history,
                                    double T) = 0;

  /// Helper to predict an elastic stress increment
  virtual Symmetric stress_increment(const Symmetric & stress,
                                     const Symmetric & D, 
                                     const Skew & W,
                                     double dt, 
                                     Lattice & lattice,
                                     const Orientation & Q,
                                     const History & history,
                                     double T) = 0;

  /// Whether this model uses the Nye tensor
  virtual bool use_nye() const;
};

/// My standard kinematic assumptions, outlined in the manual
class NEML_EXPORT StandardKinematicModel: public KinematicModel {
 public:
  /// Initialize with elastic and inelastic models
  StandardKinematicModel(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Populate a history object with the correct variables
  virtual void populate_hist(History & history) const;
  /// Initialize the history object with the starting values
  virtual void init_hist(History & history) const;

  /// Helper for external models that want a strength
  virtual double strength(const History & history, Lattice & L, double T,
                          const History & fixed) const;

  /// Anything added to this History instance will be kept fixed in the
  /// implicit integration
  virtual History decouple(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed);

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
      double T, const History & fixed) const;

  /// Helper to calculate elastic strains
  virtual Symmetric elastic_strains(const Symmetric & stress,
                                    Lattice & lattice,
                                    const Orientation & Q,
                                    const History & history,
                                    double T);
  
  /// Helper to predict an elastic stress increment
  virtual Symmetric stress_increment(const Symmetric & stress,
                                     const Symmetric & D,
                                     const Skew & W,
                                     double dt, 
                                     Lattice & lattice,
                                     const Orientation & Q,
                                     const History & history,
                                     double T);

  /// Whether this model uses the Nye tensor
  virtual bool use_nye() const;

 protected:
  std::shared_ptr<LinearElasticModel> emodel_;
  std::shared_ptr<InelasticModel> imodel_;
};

static Register<StandardKinematicModel> regStandardKinematicModel;

/// My standard kinematic assumptions with damage
class NEML_EXPORT DamagedStandardKinematicModel: public StandardKinematicModel {
 public:
  /// Initialize with elastic and inelastic models
  DamagedStandardKinematicModel(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Populate a history object with the correct variables
  virtual void populate_hist(History & history) const;
  /// Initialize the history object with the starting values
  virtual void init_hist(History & history) const;

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
      double T, const History & fixed) const;

  /// Helper to calculate elastic strains
  virtual Symmetric elastic_strains(const Symmetric & stress,
                                    Lattice & lattice,
                                    const Orientation & Q,
                                    const History & history,
                                    double T);

  /// Helper to predict an elastic stress increment
  virtual Symmetric stress_increment(const Symmetric & stress,
                                     const Symmetric & D,
                                     const Skew & W, 
                                     double dt, 
                                     Lattice & lattice,
                                     const Orientation & Q,
                                     const History & history,
                                     double T);

 protected:
  /// Names of the inelastic model internal variables
  std::vector<std::string> inames_() const;
  /// Names of the damage model internal variables
  std::vector<std::string> dnames_() const;
  /// The inelastic model internal variables
  History ihist_(const History & hist) const;
  /// The damage model internal variables
  History dhist_(const History & hist) const;
  /// Partial derivative wrt stress^\prime
  SymSymR4 d_stress_partial(
      const Symmetric & stress, const Symmetric & d,
      const Skew & w, const Orientation & Q,
      const History & history, Lattice & lattice,
      double T, const History & fixed) const;

 protected:
  std::shared_ptr<CrystalDamageModel> dmodel_;
  std::shared_ptr<AsaroInelasticity> amodel_;
};

static Register<DamagedStandardKinematicModel> regDamagedStandardKinematicModel;


} // namespace neml

#endif // KINEMATICS_H
