#pragma once

#include "singlecrystal.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include "../windows.h"

namespace neml {

class SingleCrystalModel; // forward declaration

/// Superclass of postprocessors: do something to history after update
class NEML_EXPORT CrystalPostprocessor: public NEMLObject {
 public:
  CrystalPostprocessor();

  virtual void populate_history(const Lattice & L, History & history) const = 0;
  virtual void init_history(const Lattice & L, History & history) const = 0;

  virtual void act(SingleCrystalModel & model, const Lattice &,
                   const double & T, const Symmetric & D,
                   const Skew & W, History & state,
                   const History & prev_state) = 0;
};

/// Reorients twins based on a PTR criteria
class NEML_EXPORT PTRTwinReorientation: public CrystalPostprocessor {
 public:
  PTRTwinReorientation(std::shared_ptr<Interpolate> threshold,
                       std::string prefix = "slip");

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual void populate_history(const Lattice & L, History & history) const;
  virtual void init_history(const Lattice & L, History & history) const;

  virtual void act(SingleCrystalModel & model, const Lattice &,
                   const double & T, const Symmetric & D,
                   const Skew & W, History & state,
                   const History & prev_state);

 private:
    std::shared_ptr<Interpolate> threshold_;
    std::string prefix_;
};

static Register<PTRTwinReorientation> regPTRTwinReorientation;

} // namespace neml
