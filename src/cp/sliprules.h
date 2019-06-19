#ifndef SLIPRULES_H
#define SLIPRULES_H

#include "slipharden.h"
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

class SlipHardening;

/// ABC for a slip rule
class SlipRule: public NEMLObject
{
 public:
  virtual void populate_history(History & history) const = 0;
  virtual void init_history(History & history) const = 0;

  virtual double slip(size_t g, size_t i, const Symmetric & stress, 
                      const Orientation & Q, const History & history,
                      const Lattice & L, double T) const = 0;
  virtual Symmetric d_slip_d_s(size_t g, size_t i, const Symmetric & stress, 
                               const Orientation & Q, const History & history,
                               const Lattice & L, double T) const = 0;
  virtual History
      d_slip_d_h(size_t g, size_t i, const Symmetric & stress, 
                 const Orientation & Q, const History & history,
                 const Lattice & L, double T) const = 0;

  virtual History hist_rate(const Symmetric & stress, 
                      const Orientation & Q, const History & history,
                      const Lattice & L, double T) const = 0;
  virtual History d_hist_rate_d_stress(const Symmetric & stress, 
                      const Orientation & Q, const History & history,
                      const Lattice & L, double T) const = 0;
  virtual History d_hist_rate_d_hist(const Symmetric & stress, 
                      const Orientation & Q, const History & history,
                      const Lattice & L, double T) const = 0;

};

/// All slip rules that give the system response proportional to some strength
class SlipStrengthSlipRule: public SlipRule
{
 public:
  SlipStrengthSlipRule(std::shared_ptr<SlipHardening> strength);

  virtual void populate_history(History & history) const;
  virtual void init_history(History & history) const;

  virtual double slip(size_t g, size_t i, const Symmetric & stress, 
                      const Orientation & Q, const History & history,
                      const Lattice & L, double T) const;
  virtual Symmetric d_slip_d_s(size_t g, size_t i, const Symmetric & stress, 
                               const Orientation & Q, const History & history,
                               const Lattice & L, double T) const;
  virtual History
      d_slip_d_h(size_t g, size_t i, const Symmetric & stress, 
                 const Orientation & Q, const History & history,
                 const Lattice & L, double T) const;

  virtual History hist_rate(const Symmetric & stress, 
                      const Orientation & Q, const History & history,
                      const Lattice & L, double T) const;
  virtual History d_hist_rate_d_stress(const Symmetric & stress, 
                      const Orientation & Q, const History & history,
                      const Lattice & L, double T) const;
  virtual History d_hist_rate_d_hist(const Symmetric & stress, 
                      const Orientation & Q, const History & history,
                      const Lattice & L, double T) const;

  virtual double sslip(size_t g, size_t i, double tau, double strength, 
                       double T) const = 0;
  virtual double d_sslip_dtau(size_t g, size_t i, double tau, double strength, 
                              double T) const = 0;
  virtual double d_sslip_dstrength(size_t g, size_t i, double tau, 
                                   double strength, double T) const = 0;

 private:
  std::shared_ptr<SlipHardening> strength_;
};

class PowerLawSlipRule: public SlipStrengthSlipRule
{
 public:
  PowerLawSlipRule(std::shared_ptr<SlipHardening> strength,
                   std::shared_ptr<Interpolate> gamma0, 
                   std::shared_ptr<Interpolate> n);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  virtual double sslip(size_t g, size_t i, double tau, double strength, 
                       double T) const;
  virtual double d_sslip_dtau(size_t g, size_t i, double tau, double strength, 
                              double T) const;
  virtual double d_sslip_dstrength(size_t g, size_t i, double tau, 
                                   double strength, double T) const;

 private:
  std::shared_ptr<Interpolate> gamma0_;
  std::shared_ptr<Interpolate> n_;
};

} // namespace neml

#endif // SLIPRULES_H
