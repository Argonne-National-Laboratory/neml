#ifndef INELASTICITY_H
#define INELASTICITY_H

#include "crystallography.h"
#include "sliprules.h"

#include "../objects.h"
#include "../history.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

namespace neml {

/// A inelastic model supplying D_p and W_p
//  A complete model for implicit rotations would include the derivatives wrt
//  the orientation.  Similarly it would include all the w_p derivatives.
class InelasticModel: public NEMLObject {
 public:
  virtual void populate_history(History & history) const = 0;
  virtual void init_history(History & history) const = 0;

  virtual Symmetric d_p(const Symmetric & stress, const Orientation & Q,
                        const History & history, 
                        Lattice & lattice, double T) const = 0;
  virtual SymSym d_d_p_d_stress(const Symmetric & stress, const Orientation & Q,
                                const History & history,
                                Lattice & lattice, double T) const = 0;
  virtual History d_d_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice, double T) const = 0;

  virtual History history_rate(const Symmetric & stress, const Orientation & Q,
                               const History & history,
                               Lattice & lattice, double T) const = 0;
  virtual History d_history_rate_d_stress(const Symmetric & stress, 
                                          const Orientation & Q,
                                          const History & history,
                                          Lattice & lattice,
                                          double T) const = 0;
  virtual History d_history_rate_d_history(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & lattice,
                                         double T) const = 0;

  virtual Skew w_p(const Symmetric & stress,
                   const Orientation & Q,
                   const History & history,
                   Lattice & lattice,
                   double T) const = 0;
  virtual SkewSym d_w_p_d_stress(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & lattice,
                                 double T) const = 0;
  virtual History d_w_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T) const = 0;
};

/// I don't know, maybe it might be useful
class NoInelasticity: public InelasticModel {
 public:
  NoInelasticity();
  virtual ~NoInelasticity();

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  virtual void populate_history(History & history) const;
  virtual void init_history(History & history) const;

  virtual Symmetric d_p(const Symmetric & stress,
                        const Orientation & Q,
                        const History & history,
                        Lattice & lattice,
                        double T) const;
  virtual SymSym d_d_p_d_stress(const Symmetric & stress,
                                const Orientation & Q,
                                const History & history,
                                Lattice & lattice,
                                double T) const;
  virtual History d_d_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T) const;

  virtual History history_rate(const Symmetric & stress, const Orientation & Q,
                               const History & history,
                               Lattice & lattice, double T) const;
  virtual History d_history_rate_d_stress(const Symmetric & stress, 
                                          const Orientation & Q,
                                          const History & history,
                                          Lattice & lattice,
                                          double T) const;
  virtual History d_history_rate_d_history(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & lattice,
                                         double T) const;

  virtual Skew w_p(const Symmetric & stress,
                   const Orientation & Q,
                   const History & history,
                   Lattice & lattice, double T) const;
  virtual SkewSym d_w_p_d_stress(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & lattice,
                                 double T) const;
  virtual History d_w_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T) const;
};

static Register<NoInelasticity> regNoInelasticity;

/// The classic crystal plasticity inelastic model
class AsaroInelasticity: public InelasticModel {
 public:
  AsaroInelasticity(std::shared_ptr<SlipRule> rule);
  virtual ~AsaroInelasticity();

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  virtual void populate_history(History & history) const;
  virtual void init_history(History & history) const;

  virtual Symmetric d_p(const Symmetric & stress,
                        const Orientation & Q,
                        const History & history,
                        Lattice & lattice, double T) const;
  virtual SymSym d_d_p_d_stress(const Symmetric & stress,
                                const Orientation & Q,
                                const History & history,
                                Lattice & lattice, double T) const;
  virtual History d_d_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice, double T) const;

  virtual History history_rate(const Symmetric & stress, const Orientation & Q,
                               const History & history,
                               Lattice & lattice, double T) const;
  virtual History d_history_rate_d_stress(const Symmetric & stress, 
                                          const Orientation & Q,
                                          const History & history,
                                          Lattice & lattice, double T) const;
  virtual History d_history_rate_d_history(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & lattice, double T) const;

  virtual Skew w_p(const Symmetric & stress,
                   const Orientation & Q,
                   const History & history,
                   Lattice & lattice, double T) const;
  virtual SkewSym d_w_p_d_stress(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & lattice,
                                 double T) const;
  virtual History d_w_p_d_history(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & lattice,
                                  double T) const;

 private:
  std::shared_ptr<SlipRule> rule_;
};

static Register<AsaroInelasticity> regAsaroInelasticity;

} // namespace neml

#endif // INELASTICITY_H
