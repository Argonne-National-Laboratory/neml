#pragma once

#include "slipharden.h"

#include "../objects.h"
#include "../history.h"
#include "../interpolate.h"

#include "../windows.h"

namespace neml {

/// Implementation of a single chemistry <-> size model
//  For details see Hu et al. MSE A, 2020
class NEML_EXPORT HuCocksPrecipitationModel: public HistoryNEMLObject
{
 public:
  HuCocksPrecipitationModel(ParameterSet & params);

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

  /// Actual (unscaled) f
  double f(const History & history) const;

  /// Actual (unscaled) r
  double r(const History & history) const;

  /// Actual (unscaled) N
  double N(const History & history) const;

  /// Rate vector
  std::vector<double> rate(const History & history, double T) const;

  /// Jacobian matrix
  std::vector<std::vector<double>> jac(const History & history, double T) const;

  /// The volume fraction rate
  virtual double f_rate(double f, double r, double N, double T) const;

  /// df_df
  virtual double df_df(double f, double r, double N, double T) const;

  /// df_dr
  virtual double df_dr(double f, double r, double N, double T) const;

  /// df_dN
  virtual double df_dN(double f, double r, double N, double T) const;

  /// The radius rate
  virtual double r_rate(double f, double r, double N, double T) const;

  /// dr_df
  virtual double dr_df(double f, double r, double N, double T) const;

  /// dr_dr
  virtual double dr_dr(double f, double r, double N, double T) const;

  /// dr_dN
  virtual double dr_dN(double f, double r, double N, double T) const;

  /// The number density rate
  virtual double N_rate(double f, double r, double N, double T) const;

  /// dN_df
  virtual double dN_df(double f, double r, double N, double T) const;

  /// dN_dr
  virtual double dN_dr(double f, double r, double N, double T) const;

  /// dN_dN
  virtual double dN_dN(double f, double r, double N, double T) const;

  /// Number of chemical species
  size_t nspecies() const;

  /// Concentration vector
  std::vector<double> c(double f, double T) const;

  /// Derivative of the concentration vector
  std::vector<double> dc_df(double f, double T) const;

  /// Driving force for precipitation
  double Gv(double f, double T) const;

  /// Derivative of the driving force wrt f
  double dG_df(double f, double T) const;

  /// Access the molecular volume
  double vm() const;

  /// Get the value of the volume fraction scaling term
  double fs() const {return fs_;};

  /// Get the value of the radius scaling term
  double rs() const {return rs_;};

  /// Get the value of the number density scaling term
  double Ns() const {return Ns_;};

 protected:
  /// Diffusivity
  double D_(double T) const;

  /// Check to see if we're still in the nucleation regime
  bool nucleation_(const std::vector<double> & c, double T) const;

  /// Indicator function for nucleation/ripening
  void sfn_(double f, double T, double & val, double & dval) const;

  /// The radius rate in the nucleation regime
  virtual double r_rate_nucleation_(double f, double r, double N, double T) const;

  /// dr_df in the nucleation regime
  virtual double dr_df_nucleation_(double f, double r, double N, double T) const;

  /// dr_dr in the nucleation regime
  virtual double dr_dr_nucleation_(double f, double r, double N, double T) const;

  /// dr_dN in the nucleation regime
  virtual double dr_dN_nucleation_(double f, double r, double N, double T) const;

  /// The radius rate in the ripening regime
  virtual double r_rate_ripening_(double f, double r, double N, double T) const;

  /// dr_df in the ripening regime
  virtual double dr_df_ripening_(double f, double r, double N, double T) const;

  /// dr_dr in the ripening regime
  virtual double dr_dr_ripening_(double f, double r, double N, double T) const;

  /// dr_dN in the ripening regime
  virtual double dr_dN_ripening_(double f, double r, double N, double T) const;

  /// The number density rate in the nucleation regime
  virtual double N_rate_nucleation_(double f, double r, double N, double T) const;

  /// dN_df in the nucleation regime
  virtual double dN_df_nucleation_(double f, double r, double N, double T) const;

  /// dN_dr in the nucleation regime
  virtual double dN_dr_nucleation_(double f, double r, double N, double T) const;

  /// dN_dN in the nucleation regime
  virtual double dN_dN_nucleation_(double f, double r, double N, double T) const;

  /// The number density rate in the ripening regime
  virtual double N_rate_ripening_(double f, double r, double N, double T) const;

  /// dN_df in the ripening regime
  virtual double dN_df_ripening_(double f, double r, double N, double T) const;

  /// dN_dr in the ripening regime
  virtual double dN_dr_ripening_(double f, double r, double N, double T) const;

  /// dN_dN in the ripening regime
  virtual double dN_dN_ripening_(double f, double r, double N, double T) const;

 private:
  std::vector<std::shared_ptr<Interpolate>> c0_;
  std::vector<std::shared_ptr<Interpolate>> cp_;
  std::vector<std::shared_ptr<Interpolate>> ceq_;
  double am_, N0_, Vm_, chi_, D0_, Q0_;
  std::shared_ptr<Interpolate> Cf_;
  double kboltz_, R_, Na_;
  size_t rate_;
  double f_init_, r_init_, N_init_;
  double fs_, rs_, Ns_;
  double w_;
  double vm_;
  std::vector<std::string> varnames_;
};

static Register<HuCocksPrecipitationModel> regHuCocksPrecipitationModel;

/// Standard dislocation density model, here evolving the spacing
//    See Hu and Cocks, 2020 for more details
class NEML_EXPORT DislocationSpacingHardening: public SlipHardening
{
 public:
  DislocationSpacingHardening(ParameterSet & params);

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

  /// Number of slip systems contributing
  size_t size() const;

 private:
  std::shared_ptr<Interpolate> J1_, J2_, K_;
  double L0_, a_, b_;
  std::shared_ptr<Interpolate> G_;
  std::shared_ptr<Lattice> L_;
  std::string varprefix_;
  std::vector<std::string> varnames_;
};

static Register<DislocationSpacingHardening> regDislocationSpacingHardening;

/// Full Hu and Cocks hardening model
//    See Hu and Cocks, 2020 for more details
class NEML_EXPORT HuCocksHardening: public SlipHardening
{
 public:
  HuCocksHardening(ParameterSet & params);

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
  double c_eff_(const History & history, double T) const;
  double NA_eff_(const History & history, double T) const;

 private:
  std::shared_ptr<SlipHardening> dmodel_;
  std::vector<std::shared_ptr<HuCocksPrecipitationModel>> pmodels_;
  double ap_, ac_, b_;
  std::shared_ptr<Interpolate> G_;
  std::vector<std::vector<std::string>> pnames_;
};

static Register<HuCocksHardening> regHuCocksHardening;

/// An Arrhenius slip rule ala Hu and Cocks
class NEML_EXPORT ArrheniusSlipRule: public SlipStrengthSlipRule
{
 public:
  /// Initialize with the strength object, the reference strain rate, and the
  /// rate sensitivity
  ArrheniusSlipRule(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// The slip rate definition
  virtual double scalar_sslip(size_t g, size_t i, double tau, double strength,
                              double T) const;
  /// Derivative of slip rate with respect to the resolved shear
  virtual double scalar_d_sslip_dtau(size_t g, size_t i, double tau, 
                                     double strength, double T) const;
  /// Derivative of the slip rate with respect to the strength
  virtual double scalar_d_sslip_dstrength(size_t g, size_t i, double tau,
                                          double strength, double T) const;

 private:
  double g0_, A_, B_, b_, a0_, G0_, k_;
};

static Register<ArrheniusSlipRule> regArrheniusSlipRule;


} // namespace neml
