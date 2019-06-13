#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "crystallography.h"

#include "../objects.h"
#include "../history.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

namespace neml {

class KinematicModel: public NEMLObject {
 public:
  virtual Symmetric d_p(const Symmetric & stress, const Symmetric & d,
                        const Skew & w, const History & history,
                        const Lattice & lattice) = 0;
  virtual Skew w_p(const Symmetric & stress, const Symmetric & d,
                   const Skew & w, const History & history,
                   const Lattice & lattice) = 0;
};

/// I don't know, maybe it might be useful
class NoInelasticity: public KinematicModel {
 public:
  NoInelasticity();
  virtual ~NoInelasticity();

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  virtual Symmetric d_p(const Symmetric & stress, const Symmetric & d,
                        const Skew & w, const History & history,
                        const Lattice & lattice);
  virtual Skew w_p(const Symmetric & stress, const Symmetric & d,
                   const Skew & w, const History & history,
                   const Lattice & lattice);
};

static Register<NoInelasticity> regNoInelasticity;

} // namespace neml

#endif // KINEMATICS_H
