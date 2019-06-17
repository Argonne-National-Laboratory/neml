#ifndef SLIPHARDEN_H
#define SLIPHARDEN_H

#include "sliprules.h"
#include "crystallography.h"

#include "../objects.h"
#include "../history.h"

#include "../math/tensors.h"
#include "../math/rotations.h"

#include <map>
#include <vector>
#include <string>

namespace neml {

class SlipRule; // Why would we need a forward declaration?

/// ABC for a slip hardening model
class SlipHardening: public NEMLObject
{
 public:
  /// Request whatever history you will need
  virtual void populate_history(History & history) const = 0;
  /// Setup history
  virtual void init_history(History & history) const = 0;
  
  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history) const = 0;
  /// Derivative of the map wrt to history
  virtual std::map<std::string,std::vector<double>>
      d_hist_to_tau(size_t g, size_t i, const History & history) const = 0;
  
  /// The rate of the history
  virtual History hist(const Symmetric & stress, 
                     const Orientation & Q, const History & history,
                     const Lattice & L, double T, const SlipRule & R) const = 0;
  /// Derivative of the history wrt stress
  virtual Symmetric d_hist_d_s(const Symmetric & stress, 
                               const Orientation & Q, const History & history,
                               const Lattice & L, double T,
                               const SlipRule & R) const = 0;
  /// Derivative of the history wrt the history
  virtual std::map<std::string,std::vector<double>>
      d_hist_d_h(const Symmetric & stress, 
                 const Orientation & Q,
                 const History & history,
                 const Lattice & L, 
                 double T, const SlipRule & R) const = 0;
};

/// Base for hardening models that return a single strength for all systems
class SingleSlipHardening: public SlipHardening
{
 public:
  virtual void populate_history(History & history) const;
  virtual void init_history(History & history) const;

  virtual double hist_to_tau(size_t g, size_t i, const History & history) const;
  virtual std::map<std::string,std::vector<double>>
      d_hist_to_tau(size_t g, size_t i, const History & history) const;
};

/// Base for a single hardening that depends on the sum of the slip
class SingleSSSlipHardening: public SingleSlipHardening
{
 public:
  /// The rate of the history
  virtual History hist(const Symmetric & stress, 
                     const Orientation & Q, const History & history,
                     const Lattice & L, double T, const SlipRule & R) const;
  /// Derivative of the history wrt stress
  virtual Symmetric d_hist_d_s(const Symmetric & stress, 
                               const Orientation & Q, const History & history,
                               const Lattice & L, double T,
                               const SlipRule & R) const;
  /// Derivative of the history wrt the history
  virtual std::map<std::string,std::vector<double>>
      d_hist_d_h(const Symmetric & stress, 
                 const Orientation & Q,
                 const History & history,
                 const Lattice & L, 
                 double T, const SlipRule & R) const;

  /// The scalar hardening rate
  virtual double scalar(double scalar, double T, double ss) = 0;
  /// The derivative wrt to the scalar
  virtual double d_scalar(double scalar, double T, double ss) = 0;
};

}

#endif // SLIPHARDEN_H
