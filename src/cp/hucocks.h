#pragma once

#include "slipharden.h"

#include "../objects.h"
#include "../history.h"
#include "../interpolate.h"

#include "../windows.h"

namespace neml {

/// Implementation of a single chemistry <-> size model
//  For details see Hu et al. MSE A, 2020
class NEML_EXPORT HuCocksPrecipitationModel: public NEMLObject
{
 public:
  HuCocksPrecipitationModel(std::vector<std::shared_ptr<Interpolate>> c0,
                            std::vector<std::shared_ptr<Interpolate>> cp,
                            std::vector<std::shared_ptr<Interpolate>> ceq,
                            double am, double N0, double Vm, double chi,
                            double D0, double Q0, 
                            std::shared_ptr<Interpolate> Cf,
                            double kboltz, double R, double Na,
                            size_t rate,
                            double f_init, double r_init, double N_init);

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
  virtual void populate_history(History & history) const;
  /// Setup history
  virtual void init_history(History & history) const;

  /// The volume fraction rate
  virtual double f_rate(const History & history, double T) const;

  /// df_df
  virtual double df_df(const History & history, double T) const;

  /// df_dr
  virtual double df_dr(const History & history, double T) const;

  /// df_dN
  virtual double df_dN(const History & history, double T) const;

  /// The radius rate
  virtual double r_rate(const History & history, double T) const;

  /// dr_df
  virtual double dr_df(const History & history, double T) const;

  /// dr_dr
  virtual double dr_dr(const History & history, double T) const;

  /// dr_dN
  virtual double dr_dN(const History & history, double T) const;

  /// The number density rate
  virtual double N_rate(const History & history, double T) const;

  /// dN_df
  virtual double dN_df(const History & history, double T) const;

  /// dN_dr
  virtual double dN_dr(const History & history, double T) const;

  /// dN_dN
  virtual double dN_dN(const History & history, double T) const;
  
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

 protected:
  /// Diffusivity
  double D_(double T) const;

  /// Check to see if we're still in the nucleation regime
  bool nucleation_(const std::vector<double> & c, double T) const;

 private:
  std::vector<std::shared_ptr<Interpolate>> c0_;
  std::vector<std::shared_ptr<Interpolate>> cp_;
  std::vector<std::shared_ptr<Interpolate>> ceq_;
  double am_, N0_, Vm_, chi_, D0_, Q0_;
  std::shared_ptr<Interpolate> Cf_;
  double kboltz_, R_, Na_;
  size_t rate_;
  double f_init_, r_init_, N_init_;
  double vm_;
  std::vector<std::string> varnames_;
};

static Register<HuCocksPrecipitationModel> regHuCocksPrecipitationModel;


} // namespace neml
