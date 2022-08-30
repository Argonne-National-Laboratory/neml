#ifndef SLIPHARDEN_H
#define SLIPHARDEN_H

#include "sliprules.h"
#include "crystallography.h"

#include "../objects.h"
#include "../history.h"
#include "../interpolate.h"

#include "../math/tensors.h"
#include "../math/rotations.h"
#include "../math/matrix.h"

#include "../windows.h"

#include <map>
#include <vector>
#include <string>

namespace neml {

class SlipRule; // Why would we need a forward declaration?

/// ABC for a slip hardening model
class NEML_EXPORT SlipHardening: public HistoryNEMLObject
{
 public:
  SlipHardening(ParameterSet & params);
  /// Report your variable names
  virtual std::vector<std::string> varnames() const = 0;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars) = 0;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             Lattice & L, 
                             double T, const History & fixed) const = 0;
  /// Derivative of the map wrt to history, model variables
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history, Lattice & L,
                    double T, const History & fixed) const = 0;
  /// Derivative of the map wrt to history, external variables
  virtual History
      d_hist_to_tau_ext(size_t g, size_t i, const History & history, 
                        Lattice & L, double T, const History & fixed,
                        std::vector<std::string> ext) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const = 0;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const = 0;
  /// Derivative of the history wrt the history, model variables
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const = 0;
  /// Derivative of this history wrt the history, external variables
  virtual History
      d_hist_d_h_ext(const Symmetric & stress,
                     const Orientation & Q,
                     const History & history,
                     Lattice & L,
                     double T, const SlipRule & R,
                     const History & fixed,
                     std::vector<std::string> ext) const;


  /// Whether this particular model uses the Nye tensor
  virtual bool use_nye() const;
  
  /// Helper providing a blank (zero) history
  History blank_hist() const {return cache(CacheType::BLANK);};

  enum class CacheType {BLANK, DOUBLE};

  /// Helper accessing cached objects
  History cache(CacheType type) const;

 protected:
  void init_cache_();
  std::unique_ptr<History> blank_, double_;
};

/// Fixed strength
class NEML_EXPORT FixedStrengthHardening: public SlipHardening
{
  public:
   FixedStrengthHardening(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Report your variable names
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// Request whatever history you will need
  virtual void populate_hist(History & history) const;
  /// Setup history
  virtual void init_hist(History & history) const;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             Lattice & L,
                             double T, const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history, Lattice & L,
                    double T, const History & fixed) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;

  private:
   std::vector<std::shared_ptr<Interpolate>> strengths_;
};

static Register<FixedStrengthHardening> regFixedStrengthHardening;

/// Voce model with one hardening variable per system
class NEML_EXPORT VocePerSystemHardening: public SlipHardening
{
  public:
   VocePerSystemHardening(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Report your variable names
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// Request whatever history you will need
  virtual void populate_hist(History & history) const;
  /// Setup history
  virtual void init_hist(History & history) const;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             Lattice & L,
                             double T, const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history, Lattice & L,
                    double T, const History & fixed) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;

  protected:
    size_t size_() const {return initial_.size();};
    void consistency_(Lattice & L) const;

  private:
   std::vector<double> initial_;
   std::vector<std::shared_ptr<Interpolate>> k_, sat_, m_;
   std::string varprefix_;
   std::vector<std::string> varnames_;
};

static Register<VocePerSystemHardening> regVocePerSystemHardening;

/// Frederick-Armstrong type hardening
class NEML_EXPORT FASlipHardening: public SlipHardening
{
  public:
   FASlipHardening(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Report your variable names
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// Request whatever history you will need
  virtual void populate_hist(History & history) const;
  /// Setup history
  virtual void init_hist(History & history) const;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             Lattice & L,
                             double T, const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history, Lattice & L,
                    double T, const History & fixed) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;

  protected:
    size_t size_() const {return k_.size();};
    void consistency_(Lattice & L) const;

  private:
   std::vector<std::shared_ptr<Interpolate>> k_, sat_;
   std::string varprefix_;
   std::vector<std::string> varnames_;
};

static Register<FASlipHardening> regFASlipHardening;

/// Generic linear  hardening of the form tau_i = tau_0_i + H.gamma
class NEML_EXPORT GeneralLinearHardening: public SlipHardening
{
 public:
  GeneralLinearHardening(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Report your variable names
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// Request whatever history you will need
  virtual void populate_hist(History & history) const;
  /// Setup history
  virtual void init_hist(History & history) const;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             Lattice & L,
                             double T, const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history, Lattice & L,
                    double T, const History & fixed) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;
  /// Derivative of this history wrt the history, external variables
  virtual History
      d_hist_d_h_ext(const Symmetric & stress,
                     const Orientation & Q,
                     const History & history,
                     Lattice & L,
                     double T, const SlipRule & R,
                     const History & fixed,
                     std::vector<std::string> ext) const;

 protected:
  size_t size() const {return tau_0_.size();};
  void consistency(Lattice & L) const;

 private:
  std::shared_ptr<SquareMatrix> M_;
  std::vector<double> tau_0_;
  bool absval_;
  std::string varprefix_;
  std::vector<std::string> varnames_;
};

static Register<GeneralLinearHardening> regGeneralLinearHardening;

/// Simple linear hardening, using the accumulated slip
//    \dot{h}_i = |\dot{\gamma}_i|
//    \tau_i = \tau_{0,i} + G_ij h_j
class NEML_EXPORT SimpleLinearHardening: public SlipHardening
{
 public:
  SimpleLinearHardening(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Report your variable names
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// Request whatever history you will need
  virtual void populate_hist(History & history) const;
  /// Setup history
  virtual void init_hist(History & history) const;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             Lattice & L,
                             double T, const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history, Lattice & L,
                    double T, const History & fixed) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;
  /// Derivative of this history wrt the history, external variables
  virtual History
      d_hist_d_h_ext(const Symmetric & stress,
                     const Orientation & Q,
                     const History & history,
                     Lattice & L,
                     double T, const SlipRule & R,
                     const History & fixed,
                     std::vector<std::string> ext) const;

 protected:
  size_t size() const {return tau_0_.size();};
  void consistency(Lattice & L) const;

 private:
  std::shared_ptr<SquareMatrix> G_;
  std::vector<double> tau_0_;
  std::string varprefix_;
  std::vector<std::string> varnames_;
};

static Register<SimpleLinearHardening> regSimpleLinearHardening;


class NEML_EXPORT LANLTiModel: public SlipHardening
{
 public:
  LANLTiModel(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Report your variable names
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// Request whatever history you will need
  virtual void populate_hist(History & history) const;
  /// Setup history
  virtual void init_hist(History & history) const;

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             Lattice & L,
                             double T, const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history, Lattice & L,
                    double T, const History & fixed) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;
  /// Derivative of this history wrt the history, external variables
  virtual History
      d_hist_d_h_ext(const Symmetric & stress,
                     const Orientation & Q,
                     const History & history,
                     Lattice & L,
                     double T, const SlipRule & R,
                     const History & fixed,
                     std::vector<std::string> ext) const;

 protected:
  size_t size() const {return tau_0_.size();};
  size_t nslip_() const {return k1_.size();};
  size_t ntwin_() const {return tau_0_.size()-k1_.size();};
  void consistency(Lattice & L) const;

 private:
  std::vector<std::shared_ptr<Interpolate>> tau_0_;
  std::shared_ptr<SquareMatrix> C_st_;
  std::vector<std::shared_ptr<Interpolate>> mu_, k1_, k2_;
  double X_s_;
  double inivalue_;
  std::string varprefix_, twinprefix_;
  std::vector<std::string> varnames_;
};

static Register<LANLTiModel> regLANLTiModel;


/// Slip strength rules where all systems share the same strength
class NEML_EXPORT SlipSingleHardening: public SlipHardening
{
 public:
  SlipSingleHardening(ParameterSet & params);

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g, size_t i, const History & history,
                             Lattice & L,
                             double T, const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History
      d_hist_to_tau(size_t g, size_t i, const History & history, Lattice & L,
                    double T, const History & fixed) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T, 
                          const History & fixed) const = 0;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T,
                             const History & fixed) const = 0;
};

/// Slip strength rule where all systems evolve on a single scalar strength
class NEML_EXPORT SlipSingleStrengthHardening: public SlipSingleHardening
{
 public:
  SlipSingleStrengthHardening(ParameterSet & params);

  /// Report varnames
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// Request whatever history you will need
  virtual void populate_hist(History & history) const;
  /// Setup history
  virtual void init_hist(History & history) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;

  /// Derivative of this history wrt the history, external variables
  virtual History
      d_hist_d_h_ext(const Symmetric & stress,
                     const Orientation & Q,
                     const History & history,
                     Lattice & L,
                     double T, const SlipRule & R,
                     const History & fixed,
                     std::vector<std::string> ext) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T, 
                          const History & fixed) const;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T,
                             const History & fixed) const;

  /// Set the variable's name
  void set_variable(std::string name);

  /// Static (not evolving) strength
  virtual double static_strength(double T) const = 0;

  /// Nye contribution (defaults to zero)
  double nye_contribution(const History & fixed, double T) const;

  /// Actual implementation of any Nye contribution (defaults to zero)
  virtual double nye_part(const RankTwo & nye, double T) const;

  /// Setup the scalar
  virtual double init_strength() const = 0;

  /// Scalar evolution law
  virtual double hist_rate(const Symmetric & stress, const Orientation & Q,
                           const History & history,
                           Lattice & L, double T, const SlipRule & R, 
                           const History & fixed) const = 0;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric & stress, const Orientation & Q,
                                         const History & history,
                                         Lattice & L, double T,
                                         const SlipRule & R, 
                                         const History & fixed) const = 0;
  /// Derivative of scalar law wrt the scalar
  virtual History d_hist_rate_d_hist(const Symmetric & stress, const Orientation & Q,
                                     const History & history,
                                     Lattice & L, double T,
                                     const SlipRule & R, 
                                     const History & fixed) const = 0;
  /// Derivative of the scalar law wrt all others
  virtual History d_hist_rate_d_hist_ext(const Symmetric & stress,
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & L, double T,
                                         const SlipRule & R, 
                                         const History & fixed,
                                         std::vector<std::string> ext) const = 0;
 protected:
  std::string var_name_;
};

/// Sum of individual SlipSingleStrenghHardening models (static strengths also
/// summed)
class NEML_EXPORT SumSlipSingleStrengthHardening: public SlipSingleHardening
{
 public:
  /// Initialize with a list of models
  SumSlipSingleStrengthHardening(ParameterSet & params);

  /// Report varnames
  virtual std::vector<std::string> varnames() const;
  /// Set new varnames
  virtual void set_varnames(std::vector<std::string> vars);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Request whatever history you will need
  virtual void populate_hist(History & history) const;
  /// Setup history
  virtual void init_hist(History & history) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                     const Orientation & Q, const History & history,
                     Lattice & L, double T, const SlipRule & R,
                     const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T,
                             const SlipRule & R, 
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History
      d_hist_d_h(const Symmetric & stress,
                 const Orientation & Q,
                 const History & history,
                 Lattice & L,
                 double T, const SlipRule & R,
                 const History & fixed) const;

  /// Derivative of this history wrt the history, external variables
  virtual History
      d_hist_d_h_ext(const Symmetric & stress,
                     const Orientation & Q,
                     const History & history,
                     Lattice & L,
                     double T, const SlipRule & R,
                     const History & fixed,
                     std::vector<std::string> ext) const;

  /// The scalar map
  virtual double hist_map(const History & history, double T,
                          const History & fixed) const;
  /// The derivative of the scalar map
  virtual History d_hist_map(const History & history, double T,
                             const History & fixed) const;

  /// Whether this model uses the Nye tensor
  virtual bool use_nye() const;

 private:
  size_t nmodels() const;
  const std::vector<std::shared_ptr<SlipSingleStrengthHardening>> models_;
};

static Register<SumSlipSingleStrengthHardening> regSumSlipSingleStrengthHardening;

/// Slip strength rule where the single strength evolves with sum|dg|
class NEML_EXPORT PlasticSlipHardening: public SlipSingleStrengthHardening
{
 public:
  PlasticSlipHardening(ParameterSet & params);

  /// Scalar evolution law
  virtual double hist_rate(const Symmetric & stress, const Orientation & Q,
                           const History & history, Lattice & L, double T,
                           const SlipRule & R, const History & fixed) const;
  /// Derivative of scalar law wrt stress
  virtual Symmetric d_hist_rate_d_stress(const Symmetric & stress, const Orientation & Q,
                                         const History & history, Lattice & L, double T,
                                         const SlipRule & R,
                                         const History & fixed) const;
  /// Derivative of scalar law wrt the scalar
  virtual History d_hist_rate_d_hist(const Symmetric & stress, const Orientation & Q,
                                     const History & history, Lattice & L, double T,
                                     const SlipRule & R, 
                                     const History & fixed) const;
  
  /// Derivative of the scalar law wrt all other scalars
  virtual History d_hist_rate_d_hist_ext(const Symmetric & stress, 
                                         const Orientation & Q,
                                         const History & history,
                                         Lattice & L, double T,
                                         const SlipRule & R, 
                                         const History & fixed,
                                         std::vector<std::string> ext) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice & L, double T, 
                             const History & fixed) const = 0;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice & L, double T,
                               const History & fixed) const = 0;

};

/// Everyone's favorite Voce model
class NEML_EXPORT VoceSlipHardening: public PlasticSlipHardening
{
 public:
  /// Initialize with the saturated strength, the rate constant, and a constant
  /// strength
  VoceSlipHardening(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup the scalar
  virtual double init_strength() const;

  /// Static strength
  virtual double static_strength(double T) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice & L, double T,
                             const History & fixed) const;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice & L, double T,
                               const History & fixed) const;

  /// Dynamically determine if we're going to use the Nye tensor
  virtual bool use_nye() const;

  /// Actual nye contribution
  virtual double nye_part(const RankTwo & nye, double T) const;

 private:
  std::shared_ptr<Interpolate> tau_sat_;
  std::shared_ptr<Interpolate> b_;
  std::shared_ptr<Interpolate> tau_0_;
  std::shared_ptr<Interpolate> k_;
};

static Register<VoceSlipHardening> regVoceSlipHardening;

/// The simplest of all linear hardening models
class NEML_EXPORT LinearSlipHardening: public PlasticSlipHardening
{
 public:
  /// Initialize with the saturated strength, the rate constant, and a constant
  /// strength
  LinearSlipHardening(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Setup the scalar
  virtual double init_strength() const;

  /// Static strength
  virtual double static_strength(double T) const;

  /// Prefactor
  virtual double hist_factor(double strength, Lattice & L, double T,
                             const History & fixed) const;
  /// Derivative of the prefactor
  virtual double d_hist_factor(double strength, Lattice & L, double T,
                               const History & fixed) const;

  /// Dynamically determine if we're going to use the Nye tensor
  virtual bool use_nye() const;

  /// Actual nye contribution
  virtual double nye_part(const RankTwo & nye, double T) const;

 private:
  std::shared_ptr<Interpolate> tau0_;
  std::shared_ptr<Interpolate> k1_;
  std::shared_ptr<Interpolate> k2_;
};

static Register<LinearSlipHardening> regLinearSlipHardening;

}

#endif // SLIPHARDEN_H
