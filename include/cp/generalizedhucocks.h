#pragma once

#include "slipharden.h"

#include "../history.h"
#include "../interpolate.h"
#include "../objects.h"

#include "../windows.h"

namespace neml
{

/// A chemical element in the Hu-Cocks precipitation model
class NEML_EXPORT GeneralizedHuCocksSpecies : public NEMLObject
{
public:
  GeneralizedHuCocksSpecies(ParameterSet & params);
  static std::string type() { return "GeneralizedHuCocksSpecies"; }
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params)
  {
    return neml::make_unique<GeneralizedHuCocksSpecies>(params);
  }
  static ParameterSet parameters();

  const std::string composition;
  const double c0;
  const std::shared_ptr<Interpolate> ceq;
};

static Register<GeneralizedHuCocksSpecies> regGeneralizedHuCocksSpecies;

/// A precipitate phase in the Hu-Cocks precipitation model
class NEML_EXPORT GeneralizedHuCocksPrecipitate : public HistoryNEMLObject
{
public:
  GeneralizedHuCocksPrecipitate(ParameterSet & params);
  static std::string type() { return "GeneralizedHuCocksPrecipitate"; }
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params)
  {
    return neml::make_unique<GeneralizedHuCocksPrecipitate>(params);
  }
  static ParameterSet parameters();
  virtual std::vector<std::string> varnames() const { return varnames_; }
  virtual void populate_hist(History & history) const;
  virtual void init_hist(History & history) const;
  std::map<std::string, std::shared_ptr<Interpolate>>
  species_cp_map(const std::vector<std::string> & species_names,
                 const std::vector<std::shared_ptr<Interpolate>> & cps) const;

  double f(const History & h) const { return 4. / 3. * M_PI * std::pow(r(h), 3) * N(h); }
  double d_f_d_r(const History & h) const { return 4 * M_PI * std::pow(r(h), 2) * N(h); }
  double d_f_d_N(const History & h) const { return 4. / 3. * M_PI * std::pow(r(h), 3); }
  double r(const History & history) const { return history.get<double>(varnames_[0]) * rs_; }
  double N(const History & history) const { return history.get<double>(varnames_[1]) * Ns_; }
  double rs() const { return rs_; }
  double Ns() const { return Ns_; }

  std::string composition;
  std::vector<std::shared_ptr<GeneralizedHuCocksSpecies>> species;
  std::vector<std::string> species_names;
  std::shared_ptr<GeneralizedHuCocksSpecies> rate;
  std::string rate_name;
  const std::map<std::string, std::shared_ptr<Interpolate>> cp;
  const double am, Vm, D0, Q0, N0, chi;
  const std::shared_ptr<Interpolate> Cf;

private:
  double r_init_, N_init_;
  double rs_, Ns_;
  std::vector<std::string> varnames_;
};

static Register<GeneralizedHuCocksPrecipitate> regGeneralizedHuCocksPrecipitate;

/// Implementation of the coupled chemistry <-> size model
//  For details see Hu et al. MSE A, 2020
class NEML_EXPORT GeneralizedHuCocksPrecipitationModel : public HistoryNEMLObject
{
public:
  GeneralizedHuCocksPrecipitationModel(ParameterSet & params);
  static std::string type() { return "GeneralizedHuCocksPrecipitationModel"; }
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  static ParameterSet parameters();
  virtual std::vector<std::string> varnames() const;
  virtual void populate_hist(History & history) const;
  virtual void init_hist(History & history) const;

  double diffusivity(double T,
                     const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  double concentration(double T,
                       const History & history,
                       const std::shared_ptr<GeneralizedHuCocksSpecies> & species) const;
  double
  d_concentration_d_f(double T,
                      const History & history,
                      const std::shared_ptr<GeneralizedHuCocksSpecies> & species,
                      const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  double precipitate_areal_density(const History & history) const;
  double d_precipitate_areal_density_d_r(
      const History & history,
      const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  double d_precipitate_areal_density_d_N(
      const History & history,
      const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  double
  effective_molecular_volume(const History & history,
                             const std::shared_ptr<GeneralizedHuCocksSpecies> & species) const;
  double d_effective_molecular_volume_d_f(
      const History & history,
      const std::shared_ptr<GeneralizedHuCocksSpecies> & species,
      const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  double solution_volumetric_density(double T, const History & history) const;
  double d_solution_volumetric_density_d_f(
      double T,
      const History & history,
      const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  double
  gibbs_free_energy(double T,
                    const History & history,
                    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  double d_gibbs_free_energy_d_f(
      double T,
      const History & history,
      const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
      const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const;
  std::tuple<double, double>
  growth_rate(double T,
              const History & history,
              const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  std::tuple<std::vector<double>, std::vector<double>>
  d_growth_rate(double T,
                const History & history,
                const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
                const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const;
  std::tuple<double, double>
  ripening_rate(double T,
                const History & history,
                const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  std::tuple<std::vector<double>, std::vector<double>>
  d_ripening_rate(double T,
                  const History & history,
                  const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
                  const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const;
  std::tuple<double, double>
  switching_function(double T,
                     const History & history,
                     const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
                     const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const;
  std::tuple<double, double>
  mixed_rate(double T,
             const History & history,
             const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const;
  std::tuple<std::vector<double>, std::vector<double>>
  d_mixed_rate(double T,
               const History & history,
               const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate_1,
               const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate_2) const;
  std::vector<double> rate(double T, const History & history);
  std::vector<std::vector<double>> d_rate(double T, const History & history);

  std::vector<std::shared_ptr<GeneralizedHuCocksSpecies>> species;
  std::vector<std::shared_ptr<GeneralizedHuCocksPrecipitate>> precipitates;
  const double kboltz, Na, R;
};

static Register<GeneralizedHuCocksPrecipitationModel> regGeneralizedHuCocksPrecipitationModel;

/// Full Hu and Cocks hardening model
//    See Hu and Cocks, 2020 for more details
class NEML_EXPORT GeneralizedHuCocksHardening : public SlipHardening
{
public:
  GeneralizedHuCocksHardening(ParameterSet & params);

  static std::string type() { return "GeneralizedHuCocksHardening"; }
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params)
  {
    return neml::make_unique<GeneralizedHuCocksHardening>(params);
  }
  static ParameterSet parameters();
  virtual std::vector<std::string> varnames() const;
  virtual void set_varnames(std::vector<std::string> vars) {}
  virtual void populate_hist(History & history) const;
  virtual void init_hist(History & history) const;
  const std::shared_ptr<SlipHardening> & dmodel() const { return dmodel_; }
  const std::shared_ptr<GeneralizedHuCocksPrecipitationModel> & pmodel() const { return pmodel_; }

  /// Map the set of history variables to the slip system hardening
  virtual double hist_to_tau(size_t g,
                             size_t i,
                             const History & history,
                             Lattice & L,
                             double T,
                             const History & fixed) const;
  /// Derivative of the map wrt to history
  virtual History d_hist_to_tau(size_t g,
                                size_t i,
                                const History & history,
                                Lattice & L,
                                double T,
                                const History & fixed) const;

  /// The rate of the history
  virtual History hist(const Symmetric & stress,
                       const Orientation & Q,
                       const History & history,
                       Lattice & L,
                       double T,
                       const SlipRule & R,
                       const History & fixed) const;
  /// Derivative of the history wrt stress
  virtual History d_hist_d_s(const Symmetric & stress,
                             const Orientation & Q,
                             const History & history,
                             Lattice & L,
                             double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of the history wrt the history
  virtual History d_hist_d_h(const Symmetric & stress,
                             const Orientation & Q,
                             const History & history,
                             Lattice & L,
                             double T,
                             const SlipRule & R,
                             const History & fixed) const;
  /// Derivative of this history wrt the history, external variables
  virtual History d_hist_d_h_ext(const Symmetric & stress,
                                 const Orientation & Q,
                                 const History & history,
                                 Lattice & L,
                                 double T,
                                 const SlipRule & R,
                                 const History & fixed,
                                 std::vector<std::string> ext) const;
  double ap, ac, b;
  std::shared_ptr<Interpolate> G;

private:
  std::shared_ptr<SlipHardening> dmodel_;
  std::shared_ptr<GeneralizedHuCocksPrecipitationModel> pmodel_;
};

static Register<GeneralizedHuCocksHardening> regGeneralizedHuCocksHardening;
} // namespace neml
