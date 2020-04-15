#ifndef SLIPRULES_H
#define SLIPRULES_H

#include "slipharden.h"
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

class SlipHardening;

/// Abstract base class for a slip rule
class NEML_EXPORT SlipRule: public NEMLObject
{
 public:
  /// Populate the history object with the appropriate variable names and types
  virtual void populate_history(History & history) const = 0;
  /// Actually set the history to their initial values
  virtual void init_history(History & history) const = 0;

  /// Helper for models that want an average strength
  virtual double strength(const History & history, Lattice & L, double T,
                          const History & fixed) const = 0;

  /// Slip rate on group g, system i
  virtual double slip(size_t g, size_t i, const Symmetric & stress,
                      const Orientation & Q, const History & history,
                      Lattice & L, double T, const History & fixed) const = 0;
  /// Derivative of the slip rate with respect to stress
  virtual Symmetric d_slip_d_s(size_t g, size_t i, const Symmetric & stress,
                               const Orientation & Q, const History & history,
                               Lattice & L, double T, const History & fixed) const = 0;
  /// Derivative of the slip rate with respect to history
  virtual History
      d_slip_d_h(size_t g, size_t i, const Symmetric & stress,
                 const Orientation & Q, const History & history,
                 Lattice & L, double T, const History & fixed) const = 0;

  /// History rate
  virtual History hist_rate(const Symmetric & stress,
                      const Orientation & Q, const History & history,
                      Lattice & L, double T, const History & fixed) const = 0;
  /// Derivative of the history rate with respect to the stress
  virtual History d_hist_rate_d_stress(const Symmetric & stress,
                      const Orientation & Q, const History & history,
                      Lattice & L, double T, const History & fixed) const = 0;
  /// Derivative of the history rate with respect to the history
  virtual History d_hist_rate_d_hist(const Symmetric & stress,
                      const Orientation & Q, const History & history,
                      Lattice & L, double T, const History & fixed) const = 0;

  /// Calculate the sum of the absolute value of the slip rates
  double sum_slip(const Symmetric & stress, const Orientation & Q,
                  const History & history, Lattice & L, double T,
                  const History & fixed) const;
  /// Derivative of the sum of the absolute value of the slip rates wrt stress
  Symmetric d_sum_slip_d_stress(const Symmetric & stress, const Orientation & Q,
                  const History & history, Lattice & L, double T,
                  const History & fixed) const;
  /// Derivative of the sum of the absolute value of the slip rate wrt history
  History d_sum_slip_d_hist(const Symmetric & stress, const Orientation & Q,
                  const History & history, Lattice & L, double T,
                  const History & fixed) const;

  /// Whether this model uses the Nye tensor
  virtual bool use_nye() const;
};

/// Class where all slip rules that give the system response proportional to some strength,
/// which is in turn a function of the history
class NEML_EXPORT SlipStrengthSlipRule: public SlipRule
{
 public:
  /// Initialize with the strength model
  SlipStrengthSlipRule(std::shared_ptr<SlipHardening> strength);

  /// Populate the history
  virtual void populate_history(History & history) const;
  /// Actually initialize the history
  virtual void init_history(History & history) const;

  /// Helper for models that want an average strength
  virtual double strength(const History & history, Lattice & L, double T, 
                          const History & fixed) const;

  /// Slip rate on group g, system i
  virtual double slip(size_t g, size_t i, const Symmetric & stress,
                      const Orientation & Q, const History & history,
                      Lattice & L, double T,
                      const History & fixed) const;
  /// Derivative of slip rate with respect to stress
  virtual Symmetric d_slip_d_s(size_t g, size_t i, const Symmetric & stress,
                               const Orientation & Q, const History & history,
                               Lattice & L, double T,
                               const History & fixed) const;
  /// Derivative of slip rate with respect to history
  virtual History
      d_slip_d_h(size_t g, size_t i, const Symmetric & stress,
                 const Orientation & Q, const History & history,
                 Lattice & L, double T, const History & fixed) const;

  /// History evolution equations
  virtual History hist_rate(const Symmetric & stress,
                      const Orientation & Q, const History & history,
                      Lattice & L, double T, const History & fixed) const;
  /// Derivative of the history rate with respect to stress
  virtual History d_hist_rate_d_stress(const Symmetric & stress,
                      const Orientation & Q, const History & history,
                      Lattice & L, double T, const History & fixed) const;
  /// Derivative of the history rate with respect to the history
  virtual History d_hist_rate_d_hist(const Symmetric & stress,
                      const Orientation & Q, const History & history,
                      Lattice & L, double T, const History & fixed) const;

  /// The slip rate on group g, system i given the resolved shear, the strength,
  /// and temperature
  virtual double sslip(size_t g, size_t i, double tau, double strength,
                       double T) const = 0;
  /// Derivative of slip rate with respect to the resolved shear
  virtual double d_sslip_dtau(size_t g, size_t i, double tau, double strength,
                              double T) const = 0;
  /// Derivative of the slip rate with respect to the strength
  virtual double d_sslip_dstrength(size_t g, size_t i, double tau,
                                   double strength, double T) const = 0;

  /// Whether this model uses the Nye tensor
  virtual bool use_nye() const;

 private:
  std::shared_ptr<SlipHardening> strength_;
};

/// The standard power law slip strength/rate relation
class NEML_EXPORT PowerLawSlipRule: public SlipStrengthSlipRule
{
 public:
  /// Initialize with the strength object, the reference strain rate, and the
  /// rate sensitivity
  PowerLawSlipRule(std::shared_ptr<SlipHardening> strength,
                   std::shared_ptr<Interpolate> gamma0,
                   std::shared_ptr<Interpolate> n);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// The slip rate definition
  virtual double sslip(size_t g, size_t i, double tau, double strength,
                       double T) const;
  /// Derivative of slip rate with respect to the resolved shear
  virtual double d_sslip_dtau(size_t g, size_t i, double tau, double strength,
                              double T) const;
  /// Derivative of the slip rate with respect to the strength
  virtual double d_sslip_dstrength(size_t g, size_t i, double tau,
                                   double strength, double T) const;

 private:
  std::shared_ptr<Interpolate> gamma0_;
  std::shared_ptr<Interpolate> n_;
};

static Register<PowerLawSlipRule> regPowerLawSlipRule;

} // namespace neml

#endif // SLIPRULES_H
