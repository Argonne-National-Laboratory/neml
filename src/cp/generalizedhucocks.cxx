#include "cp/generalizedhucocks.h"
#include "parse.h"

namespace neml
{

GeneralizedHuCocksSpecies::GeneralizedHuCocksSpecies(ParameterSet & params)
  : NEMLObject(params),
    composition(params.get_parameter<std::string>("composition")),
    c0(params.get_parameter<double>("c0")),
    ceq(params.get_object_parameter<Interpolate>("ceq"))
{
}

ParameterSet
GeneralizedHuCocksSpecies::parameters()
{
  ParameterSet pset(GeneralizedHuCocksSpecies::type());

  pset.add_parameter<std::string>("composition");
  pset.add_parameter<double>("c0");
  pset.add_parameter<NEMLObject>("ceq");

  return pset;
}

GeneralizedHuCocksPrecipitate::GeneralizedHuCocksPrecipitate(ParameterSet & params)
  : HistoryNEMLObject(params),
    // public
    composition(params.get_parameter<std::string>("composition")),
    species_names(split_string(params.get_parameter<std::string>("species"))),
    rate_name(params.get_parameter<std::string>("rate")),
    cp(species_cp_map(species_names, params.get_object_parameter_vector<Interpolate>("cp"))),
    am(params.get_parameter<double>("am")),
    Vm(params.get_parameter<double>("Vm")),
    D0(params.get_parameter<double>("D0")),
    Q0(params.get_parameter<double>("Q0")),
    N0(params.get_parameter<double>("N0")),
    chi(params.get_parameter<double>("chi")),
    Cf(params.get_object_parameter<Interpolate>("Cf")),
    // private
    r_init_(params.get_parameter<double>("r_init")),
    N_init_(params.get_parameter<double>("N_init")),
    rs_(params.get_parameter<double>("rs")),
    Ns_(params.get_parameter<double>("Ns")),
    varnames_({composition + "_r", composition + "_N"})
{
}

ParameterSet
GeneralizedHuCocksPrecipitate::parameters()
{
  ParameterSet pset(GeneralizedHuCocksPrecipitate::type());

  pset.add_parameter<std::string>("composition");
  pset.add_parameter<std::vector<NEMLObject>>("cp");
  pset.add_parameter<double>("am");
  pset.add_parameter<double>("Vm");
  pset.add_parameter<double>("D0");
  pset.add_parameter<double>("Q0");
  pset.add_parameter<double>("N0");
  pset.add_parameter<double>("chi");
  pset.add_parameter<NEMLObject>("Cf");
  pset.add_parameter<std::string>("species");
  pset.add_parameter<std::string>("rate");

  pset.add_optional_parameter<double>("r_init", 1.0e-9);
  pset.add_optional_parameter<double>("N_init", 1.0e11);
  pset.add_optional_parameter<double>("rs", 1.0e-9);
  pset.add_optional_parameter<double>("Ns", 1.0e12);

  return pset;
}

void
GeneralizedHuCocksPrecipitate::populate_hist(History & history) const
{
  for (auto vn : varnames_)
    history.add<double>(vn);
}

void
GeneralizedHuCocksPrecipitate::init_hist(History & history) const
{
  history.get<double>(varnames_[0]) = r_init_ / rs_;
  history.get<double>(varnames_[1]) = N_init_ / Ns_;
}

std::map<std::string, std::shared_ptr<Interpolate>>
GeneralizedHuCocksPrecipitate::species_cp_map(
    const std::vector<std::string> & species_names,
    const std::vector<std::shared_ptr<Interpolate>> & cps) const
{
  std::map<std::string, std::shared_ptr<Interpolate>> res;
  std::transform(species_names.begin(),
                 species_names.end(),
                 cps.begin(),
                 std::inserter(res, res.end()),
                 [](std::string k, std::shared_ptr<Interpolate> v)
                 { return std::make_pair(k, v); });
  return res;
}

GeneralizedHuCocksPrecipitationModel::GeneralizedHuCocksPrecipitationModel(ParameterSet & params)
  : HistoryNEMLObject(params),
    // public
    species(params.get_object_parameter_vector<GeneralizedHuCocksSpecies>("species")),
    precipitates(params.get_object_parameter_vector<GeneralizedHuCocksPrecipitate>("precipitates")),
    kboltz(params.get_parameter<double>("kboltz")),
    Na(params.get_parameter<double>("Na")),
    R(params.get_parameter<double>("R"))
{
}

std::unique_ptr<NEMLObject>
GeneralizedHuCocksPrecipitationModel::initialize(ParameterSet & params)
{
  auto model = neml::make_unique<GeneralizedHuCocksPrecipitationModel>(params);

  for (auto & precipitate : model->precipitates)
  {
    // Let each precipitate know where to find the species its composed of.
    // We can only do this here as one species may go into multiple precipitates.
    for (auto & species_name : precipitate->species_names)
    {
      bool found = false;
      for (auto & species : model->species)
        if (species->composition == species_name)
        {
          precipitate->species.push_back(species);
          found = true;
          break;
        }
      if (!found)
        throw std::runtime_error("Unrecognized species compositon " + species_name);
    }
    // Also let each precipitate know the pointer to its rate limiting species.
    // We can only do this here as one species may go into multiple precipitates.

    bool found = false;
    for (auto & species : model->species)
      if (species->composition == precipitate->rate_name)
      {
        precipitate->rate = species;
        found = true;
        break;
      }
    if (!found)
      throw std::runtime_error("Unrecognized rate-limiting species compositon " +
                               precipitate->rate_name);
  }

  return model;
}

ParameterSet
GeneralizedHuCocksPrecipitationModel::parameters()
{
  ParameterSet pset(GeneralizedHuCocksPrecipitationModel::type());

  pset.add_parameter<std::vector<NEMLObject>>("species");
  pset.add_parameter<std::vector<NEMLObject>>("precipitates");
  pset.add_optional_parameter<double>("kboltz", 1.3806485e-23);
  pset.add_optional_parameter<double>("Na", 6.02e23);
  pset.add_optional_parameter<double>("R", 8.31462);

  return pset;
}

std::vector<std::string>
GeneralizedHuCocksPrecipitationModel::varnames() const
{
  std::vector<std::string> names;
  for (auto & p : precipitates)
  {
    const auto & pn = p->varnames();
    names.insert(names.end(), pn.begin(), pn.end());
  }
  return names;
}

void
GeneralizedHuCocksPrecipitationModel::populate_hist(History & history) const
{
  for (auto p : precipitates)
    p->populate_hist(history);
}

void
GeneralizedHuCocksPrecipitationModel::init_hist(History & history) const
{
  for (auto p : precipitates)
    p->init_hist(history);
}

double
GeneralizedHuCocksPrecipitationModel::diffusivity(
    double T, const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  return precipitate->D0 * std::exp(-precipitate->Q0 / (R * T));
}

double
GeneralizedHuCocksPrecipitationModel::concentration(
    double T, const History & history, const std::shared_ptr<GeneralizedHuCocksSpecies> & ss) const
{
  double c_total = ss->c0;
  double f_total = 1;
  for (auto & p : precipitates)
    for (auto & s : p->species)
      if (s == ss)
      {
        c_total -= p->f(history) * p->cp.at(s->composition)->value(T);
        f_total -= p->f(history);
      }
  double c = c_total / f_total;
  return c > 1e-12 ? c : 1e-12;
}

double
GeneralizedHuCocksPrecipitationModel::d_concentration_d_f(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksSpecies> & ss,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & pp) const
{
  if (pp->cp.count(ss->composition) == 0)
    return 0;

  double d_c_total = -pp->cp.at(ss->composition)->value(T);
  double d_f_total = -1;

  double c_total = ss->c0;
  double f_total = 1;
  for (auto & p : precipitates)
    for (auto & s : p->species)
      if (s->composition == ss->composition)
      {
        c_total -= p->f(history) * p->cp.at(s->composition)->value(T);
        f_total -= p->f(history);
      }

  double c = c_total / f_total;
  double d_c_d_f = c > 1e-12 ? d_c_total / f_total - c_total / f_total / f_total * d_f_total : 0;
  return d_c_d_f;
}

double
GeneralizedHuCocksPrecipitationModel::precipitate_areal_density(const History & history) const
{
  double res = 0;
  for (const auto & p : precipitates)
    res += 2 * p->r(history) * p->N(history);
  return res;
}

double
GeneralizedHuCocksPrecipitationModel::d_precipitate_areal_density_d_r(
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  return 2 * precipitate->N(history);
}

double
GeneralizedHuCocksPrecipitationModel::d_precipitate_areal_density_d_N(
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  return 2 * precipitate->r(history);
}

double
GeneralizedHuCocksPrecipitationModel::effective_molecular_volume(
    const History & history, const std::shared_ptr<GeneralizedHuCocksSpecies> & species) const
{
  double fv_total = 0;
  double f_total = 0;
  for (const auto & p : precipitates)
    for (const auto & s : p->species)
      if (s == species)
      {
        fv_total += p->f(history) * p->Vm / Na;
        f_total += p->f(history);
      }
  return fv_total / f_total;
}

double
GeneralizedHuCocksPrecipitationModel::d_effective_molecular_volume_d_f(
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksSpecies> & species,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  double fv_total = 0;
  double f_total = 0;
  for (const auto & p : precipitates)
    for (const auto & s : p->species)
      if (s->composition == species->composition)
      {
        fv_total += p->f(history) * p->Vm / Na;
        f_total += p->f(history);
      }

  double d_fv_total_d_f = 0;
  double d_f_total_d_f = 0;
  for (const auto & s : precipitate->species)
    if (s->composition == species->composition)
    {
      d_fv_total_d_f = precipitate->Vm / Na;
      d_f_total_d_f = 1;
      break;
    }

  return d_fv_total_d_f / f_total - fv_total / f_total / f_total * d_f_total_d_f;
}

double
GeneralizedHuCocksPrecipitationModel::solution_volumetric_density(double T,
                                                                  const History & history) const
{
  double res = 0;
  for (const auto & s : species)
    res += concentration(T, history, s) / effective_molecular_volume(history, s);
  return res;
}

double
GeneralizedHuCocksPrecipitationModel::d_solution_volumetric_density_d_f(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  double res = 0;
  for (const auto & s : species)
  {
    double c = concentration(T, history, s);
    double d_c_d_f = d_concentration_d_f(T, history, s, precipitate);
    double vm = effective_molecular_volume(history, s);
    double d_vm_d_f = d_effective_molecular_volume_d_f(history, s, precipitate);
    res += d_c_d_f / vm - c / vm / vm * d_vm_d_f;
  }
  return res;
}

double
GeneralizedHuCocksPrecipitationModel::gibbs_free_energy(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  double vm = precipitate->Vm / Na;
  double Gv = 0;
  for (auto & s : precipitate->species)
    Gv += -std::log(concentration(T, history, s) / s->ceq->value(T));
  Gv *= kboltz * T / vm;
  return Gv;
}

double
GeneralizedHuCocksPrecipitationModel::d_gibbs_free_energy_d_f(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const
{
  double vm = precipitate1->Vm / Na;
  double d_Gv_d_f = 0;
  for (auto & s : precipitate1->species)
    d_Gv_d_f +=
        -1 / concentration(T, history, s) * d_concentration_d_f(T, history, s, precipitate2);
  d_Gv_d_f *= kboltz * T / vm;
  return d_Gv_d_f;
}

std::tuple<double, double>
GeneralizedHuCocksPrecipitationModel::growth_rate(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  double r = precipitate->r(history);
  double N = precipitate->N(history);
  double N0 = precipitate->N0;
  double vm = precipitate->Vm / Na;
  double am = precipitate->am;
  double chi = precipitate->chi;
  auto rate = precipitate->rate;
  double ceq = rate->ceq->value(T);
  double cp = precipitate->cp.at(rate->composition)->value(T);

  double D = diffusivity(T, precipitate);
  double c = concentration(T, history, rate);
  double Zbeta = 2 * vm * D * c / std::pow(am, 4) * std::sqrt(chi / kboltz / T);
  double Gv = gibbs_free_energy(T, history, precipitate);
  double Gs = 16 * M_PI * std::pow(chi, 3) / 3 / std::pow(Gv, 2);
  double rc = -2 * chi / Gv;

  double N_dot = N0 * Zbeta * std::exp(-Gs / kboltz / T);
  double r_dot = D / r * (c - ceq) / (cp - ceq) + N_dot / N * (rc - r);

  return {r_dot, N_dot};
}

std::tuple<std::vector<double>, std::vector<double>>
GeneralizedHuCocksPrecipitationModel::d_growth_rate(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const
{
  double r = precipitate1->r(history);
  double N = precipitate1->N(history);
  double N0 = precipitate1->N0;
  double vm = precipitate1->Vm / Na;
  double am = precipitate1->am;
  double chi = precipitate1->chi;
  auto rate = precipitate1->rate;
  double ceq = rate->ceq->value(T);
  double cp = precipitate1->cp.at(rate->composition)->value(T);

  double D = diffusivity(T, precipitate1);
  double c = concentration(T, history, rate);
  double Zbeta = 2 * vm * D * c / std::pow(am, 4) * std::sqrt(chi / kboltz / T);
  double Gv = gibbs_free_energy(T, history, precipitate1);
  double Gs = 16 * M_PI * std::pow(chi, 3) / 3 / std::pow(Gv, 2);
  double rc = -2 * chi / Gv;

  double N_dot = N0 * Zbeta * std::exp(-Gs / kboltz / T);

  double d_N_dot_d_Zbeta = N0 * std::exp(-Gs / kboltz / T);
  double d_Zbeta_d_c = 2 * vm * D / std::pow(am, 4) * std::sqrt(chi / kboltz / T);
  double d_c_d_f = d_concentration_d_f(T, history, rate, precipitate2);
  double d_N_dot_d_Gs = N0 * Zbeta * std::exp(-Gs / kboltz / T) * -1 / kboltz / T;
  double d_Gs_d_Gv = -32 * M_PI * std::pow(chi, 3) / 3 / std::pow(Gv, 3);
  double d_Gv_d_f = d_gibbs_free_energy_d_f(T, history, precipitate1, precipitate2);
  double d_N_dot_d_f =
      d_N_dot_d_Zbeta * d_Zbeta_d_c * d_c_d_f + d_N_dot_d_Gs * d_Gs_d_Gv * d_Gv_d_f;
  double d_N_dot_d_r = d_N_dot_d_f * precipitate2->d_f_d_r(history);
  double d_N_dot_d_N = d_N_dot_d_f * precipitate2->d_f_d_N(history);

  double d_r_dot_d_rc = N_dot / N;
  double d_rc_d_Gv = 2 * chi / Gv / Gv;
  double d_r_dot_d_c = D / r / (cp - ceq);
  double d_r_dot_d_f =
      d_r_dot_d_c * d_c_d_f + d_r_dot_d_rc * d_rc_d_Gv * d_Gv_d_f + d_N_dot_d_f / N * (rc - r);
  double d_r_dot_d_r = d_r_dot_d_f * precipitate2->d_f_d_r(history);
  double d_r_dot_d_N = d_r_dot_d_f * precipitate2->d_f_d_N(history);
  d_r_dot_d_r += precipitate1 == precipitate2 ? -D / r / r * (c - ceq) / (cp - ceq) - N_dot / N : 0;
  d_r_dot_d_N += precipitate1 == precipitate2 ? -N_dot / N / N * (rc - r) : 0;

  return {{d_r_dot_d_r, d_r_dot_d_N}, {d_N_dot_d_r, d_N_dot_d_N}};
}

std::tuple<double, double>
GeneralizedHuCocksPrecipitationModel::ripening_rate(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  double r = precipitate->r(history);
  double N = precipitate->N(history);
  double Cf = precipitate->Cf->value(T);
  double chi = precipitate->chi;
  double Vm = precipitate->Vm;
  auto rate = precipitate->rate;

  double D = diffusivity(T, precipitate);
  double c = concentration(T, history, rate);
  double M = Cf * 8 * chi * Vm * D * c / 9 / R / T;

  double r_dot = M / 3 / r / r;
  double N_dot = -3 * N / r * r_dot;

  return {r_dot, N_dot};
}

std::tuple<std::vector<double>, std::vector<double>>
GeneralizedHuCocksPrecipitationModel::d_ripening_rate(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const
{
  double r = precipitate1->r(history);
  double N = precipitate1->N(history);
  double Cf = precipitate1->Cf->value(T);
  double chi = precipitate1->chi;
  double Vm = precipitate1->Vm;
  auto rate = precipitate1->rate;

  double D = diffusivity(T, precipitate1);
  double c = concentration(T, history, rate);
  double M = Cf * 8 * chi * Vm * D * c / 9 / R / T;

  double r_dot = M / 3 / r / r;

  double d_r_dot_d_M = 1. / 3. / r / r;
  double d_M_d_c = Cf * 8 * chi * Vm * D / 9 / R / T;
  double d_c_d_f = d_concentration_d_f(T, history, rate, precipitate2);
  double d_r_dot_d_f = d_r_dot_d_M * d_M_d_c * d_c_d_f;
  double d_r_dot_d_r = d_r_dot_d_f * precipitate2->d_f_d_r(history);
  double d_r_dot_d_N = d_r_dot_d_f * precipitate2->d_f_d_N(history);
  d_r_dot_d_r += precipitate1 == precipitate2 ? -2 * M / 3 / r / r / r : 0;

  double d_N_dot_d_r_dot = -3 * N / r;
  double d_N_dot_d_f = d_N_dot_d_r_dot * d_r_dot_d_f;
  double d_N_dot_d_r = d_N_dot_d_f * precipitate2->d_f_d_r(history);
  double d_N_dot_d_N = d_N_dot_d_f * precipitate2->d_f_d_N(history);
  d_N_dot_d_r +=
      (precipitate1 == precipitate2 ? 3 * N / r / r * r_dot : 0) + d_N_dot_d_r_dot * d_r_dot_d_r;
  d_N_dot_d_N += precipitate1 == precipitate2 ? -3 / r * r_dot : 0;

  return {{d_r_dot_d_r, d_r_dot_d_N}, {d_N_dot_d_r, d_N_dot_d_N}};
}

std::tuple<double, double>
GeneralizedHuCocksPrecipitationModel::switching_function(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const
{
  double fi = 0;
  double d_fi_d_f = 0;
  for (auto & s : precipitate1->species)
  {
    double c = concentration(T, history, s);
    double fi_s = (c - s->c0) / (s->ceq->value(T) - s->c0);
    if (fi_s > fi)
    {
      fi = fi_s;
      d_fi_d_f = d_concentration_d_f(T, history, s, precipitate2) / (s->ceq->value(T) - s->c0);
    }
  }
  if (fi > 1)
  {
    fi = 1;
    d_fi_d_f = 0;
  }
  return {fi, d_fi_d_f};
}

std::tuple<double, double>
GeneralizedHuCocksPrecipitationModel::mixed_rate(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate) const
{
  const auto [fi, d_fi_d_f] = switching_function(T, history, precipitate, precipitate);

  auto [r_dot_growth, N_dot_growth] = growth_rate(T, history, precipitate);
  auto [r_dot_ripening, N_dot_ripening] = ripening_rate(T, history, precipitate);
  double r_dot = (1 - fi) * r_dot_growth + fi * r_dot_ripening;
  double N_dot = (1 - fi) * N_dot_growth + fi * N_dot_ripening;

  return {r_dot, N_dot};
}

std::tuple<std::vector<double>, std::vector<double>>
GeneralizedHuCocksPrecipitationModel::d_mixed_rate(
    double T,
    const History & history,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate1,
    const std::shared_ptr<GeneralizedHuCocksPrecipitate> & precipitate2) const
{
  const auto [fi, d_fi_d_f] = switching_function(T, history, precipitate1, precipitate2);

  auto [r_dot_growth, N_dot_growth] = growth_rate(T, history, precipitate1);
  auto [r_dot_ripening, N_dot_ripening] = ripening_rate(T, history, precipitate1);

  auto [d_r_dot_growth, d_N_dot_growth] = d_growth_rate(T, history, precipitate1, precipitate2);
  auto [d_r_dot_ripening, d_N_dot_ripening] =
      d_ripening_rate(T, history, precipitate1, precipitate2);

  double d_r_dot_d_f = d_fi_d_f * (r_dot_ripening - r_dot_growth);
  double d_r_dot_d_r = (1 - fi) * d_r_dot_growth[0] + fi * d_r_dot_ripening[0] +
                       d_r_dot_d_f * precipitate2->d_f_d_r(history);
  double d_r_dot_d_N = (1 - fi) * d_r_dot_growth[1] + fi * d_r_dot_ripening[1] +
                       d_r_dot_d_f * precipitate2->d_f_d_N(history);

  double d_N_dot_d_f = d_fi_d_f * (N_dot_ripening - N_dot_growth);
  double d_N_dot_d_r = (1 - fi) * d_N_dot_growth[0] + fi * d_N_dot_ripening[0] +
                       d_N_dot_d_f * precipitate2->d_f_d_r(history);
  double d_N_dot_d_N = (1 - fi) * d_N_dot_growth[1] + fi * d_N_dot_ripening[1] +
                       d_N_dot_d_f * precipitate2->d_f_d_N(history);

  return {{d_r_dot_d_r, d_r_dot_d_N}, {d_N_dot_d_r, d_N_dot_d_N}};
}

std::vector<double>
GeneralizedHuCocksPrecipitationModel::rate(double T, const History & history)
{
  std::vector<double> res(precipitates.size() * 2);
  for (unsigned int i = 0; i < precipitates.size(); i++)
  {
    const auto & pi = precipitates[i];
    const auto [r_dot, N_dot] = mixed_rate(T, history, pi);
    res[2 * i] = r_dot / pi->rs();
    res[2 * i + 1] = N_dot / pi->Ns();
  }
  return res;
}

std::vector<std::vector<double>>
GeneralizedHuCocksPrecipitationModel::d_rate(double T, const History & history)
{
  std::vector<std::vector<double>> jac(precipitates.size() * 2);
  for (unsigned int i = 0; i < precipitates.size(); i++)
  {
    const auto & pi = precipitates[i];
    jac[2 * i].resize(precipitates.size() * 2);
    jac[2 * i + 1].resize(precipitates.size() * 2);
    for (unsigned int j = 0; j < precipitates.size(); j++)
    {
      const auto & pj = precipitates[j];
      const auto [d_r_dot, d_N_dot] = d_mixed_rate(T, history, pi, pj);
      jac[2 * i][2 * j] = d_r_dot[0] / pi->rs() * pj->rs();
      jac[2 * i][2 * j + 1] = d_r_dot[1] / pi->rs() * pj->Ns();

      jac[2 * i + 1][2 * j] = d_N_dot[0] / pi->Ns() * pj->rs();
      jac[2 * i + 1][2 * j + 1] = d_N_dot[1] / pi->Ns() * pj->Ns();
    }
  }
  return jac;
}

GeneralizedHuCocksHardening::GeneralizedHuCocksHardening(ParameterSet & params)
  : SlipHardening(params),
    ap(params.get_parameter<double>("ap")),
    ac(params.get_parameter<double>("ac")),
    b(params.get_parameter<double>("b")),
    G(params.get_object_parameter<Interpolate>("G")),
    dmodel_(params.get_object_parameter<SlipHardening>("dmodel")),
    pmodel_(params.get_object_parameter<GeneralizedHuCocksPrecipitationModel>("pmodel"))
{
  init_cache_();
}

ParameterSet
GeneralizedHuCocksHardening::parameters()
{
  ParameterSet pset(GeneralizedHuCocksHardening::type());

  pset.add_parameter<NEMLObject>("dmodel");
  pset.add_parameter<NEMLObject>("pmodel");
  pset.add_parameter<double>("ap");
  pset.add_parameter<double>("ac");
  pset.add_parameter<double>("b");
  pset.add_parameter<NEMLObject>("G");

  return pset;
}

std::vector<std::string>
GeneralizedHuCocksHardening::varnames() const
{
  auto names = dmodel_->varnames();
  auto pnames = pmodel_->varnames();
  names.insert(names.end(), pnames.begin(), pnames.end());
  return names;
}

void
GeneralizedHuCocksHardening::populate_hist(History & history) const
{
  dmodel_->populate_hist(history);
  pmodel_->populate_hist(history);
}

void
GeneralizedHuCocksHardening::init_hist(History & history) const
{
  history.zero();
  dmodel_->init_hist(history);
  pmodel_->init_hist(history);
}

double
GeneralizedHuCocksHardening::hist_to_tau(
    size_t g, size_t i, const History & history, Lattice & L, double T, const History & fixed) const
{
  double svd = pmodel_->solution_volumetric_density(T, history);
  double pad = pmodel_->precipitate_areal_density(history);

  double tau_d = dmodel_->hist_to_tau(g, i, history, L, T, fixed);
  double tau_p = ap * G->value(T) * b * std::sqrt(pad);
  double tau_c = ac * G->value(T) * b * std::sqrt(svd * b);

  return std::sqrt(tau_d * tau_d + tau_p * tau_p) + tau_c;
}

History
GeneralizedHuCocksHardening::d_hist_to_tau(
    size_t g, size_t i, const History & history, Lattice & L, double T, const History & fixed) const
{
  History res = cache(CacheType::DOUBLE).zero();

  // Commonly-used things
  double svd = pmodel_->solution_volumetric_density(T, history);
  double pad = pmodel_->precipitate_areal_density(history);

  double tau_d = dmodel_->hist_to_tau(g, i, history, L, T, fixed);
  double tau_p = ap * G->value(T) * b * std::sqrt(pad);
  double tau_c = ac * G->value(T) * b * std::sqrt(svd * b);
  double tau_s = std::sqrt(tau_p * tau_p + tau_d * tau_d);

  // First block: dislocation terms
  History dd = dmodel_->d_hist_to_tau(g, i, history, L, T, fixed);
  dd.scalar_multiply(tau_d / tau_s);
  std::copy(dd.rawptr(), dd.rawptr() + dd.size(), res.rawptr());

  // For each precipitation model
  for (const auto & p : pmodel_->precipitates)
  {
    double d_svd_d_f = pmodel_->d_solution_volumetric_density_d_f(T, history, p);
    double d_pad_d_r = pmodel_->d_precipitate_areal_density_d_r(history, p);
    double d_pad_d_N = pmodel_->d_precipitate_areal_density_d_N(history, p);

    // Third block: r
    res.get<double>(p->varnames()[0]) = (0.5 * tau_c / svd * d_svd_d_f * p->d_f_d_r(history) +
                                         0.5 * tau_p * tau_p / tau_s / pad * d_pad_d_r) *
                                        p->rs();

    // Fourth block: N
    res.get<double>(p->varnames()[1]) = (0.5 * tau_c / svd * d_svd_d_f * p->d_f_d_N(history) +
                                         0.5 * tau_p * tau_p / tau_s / pad * d_pad_d_N) *
                                        p->Ns();
  }

  return res;
}

History
GeneralizedHuCocksHardening::hist(const Symmetric & stress,
                                  const Orientation & Q,
                                  const History & history,
                                  Lattice & L,
                                  double T,
                                  const SlipRule & R,
                                  const History & fixed) const
{
  // Vector of results
  History res = blank_hist().zero();

  // History is easy, we just concatenate the dislocation model hardening
  // followed by each precipitation hardening
  auto h1 = dmodel_->hist(stress, Q, history, L, T, R, fixed);
  std::copy(h1.rawptr(), h1.rawptr() + h1.size(), res.start_loc(dmodel_->varnames()[0]));

  auto rate = pmodel_->rate(T, history);
  for (unsigned int i = 0; i < pmodel_->varnames().size(); i++)
    res.get<double>(pmodel_->varnames()[i]) = rate[i];

  return res;
}

History
GeneralizedHuCocksHardening::d_hist_d_s(const Symmetric & stress,
                                        const Orientation & Q,
                                        const History & history,
                                        Lattice & L,
                                        double T,
                                        const SlipRule & R,
                                        const History & fixed) const
{
  // This could be non-zero
  History res = dmodel_->d_hist_d_s(stress, Q, history, L, T, R, fixed);

  // These are zero
  History pm;
  pmodel_->populate_hist(pm);
  res.add_union(pm.derivative<Symmetric>().zero());

  return res;
}

History
GeneralizedHuCocksHardening::d_hist_d_h(const Symmetric & stress,
                                        const Orientation & Q,
                                        const History & history,
                                        Lattice & L,
                                        double T,
                                        const SlipRule & R,
                                        const History & fixed) const
{
  // Start with the self-derivative
  auto res = dmodel_->d_hist_d_h(stress, Q, history, L, T, R, fixed);

  // pmodel/dmodel cross term
  History phist;
  pmodel_->populate_hist(phist);
  res.add_union(phist.history_derivative(dmodel_->blank_hist()).zero());

  // dmodel/pmodel cross term
  // res.add_union(dmodel_->blank_hist().history_derivative(phists[i]).zero());

  // pmodel jacobian
  auto jac = pmodel_->d_rate(T, history);
  for (size_t i = 0; i < pmodel_->varnames().size(); i++)
    for (size_t j = 0; j < pmodel_->varnames().size(); j++)
    {
      res.add<double>(pmodel_->varnames()[i] + "_" + pmodel_->varnames()[j]);
      res.get<double>(pmodel_->varnames()[i] + "_" + pmodel_->varnames()[j]) = jac[i][j];
    }

  // Reorder...
  std::vector<std::string> order;
  for (auto n1 : varnames())
    for (auto n2 : varnames())
      order.push_back(n1 + "_" + n2);
  res.reorder(order);

  return res;
}

History
GeneralizedHuCocksHardening::d_hist_d_h_ext(const Symmetric & stress,
                                            const Orientation & Q,
                                            const History & history,
                                            Lattice & L,
                                            double T,
                                            const SlipRule & R,
                                            const History & fixed,
                                            std::vector<std::string> ext) const
{
  History res = blank_hist().history_derivative(history.subset(ext)).zero();

  // The only potential non-zero is from the dmodel
  History h1 = dmodel_->d_hist_d_h_ext(stress, Q, history, L, T, R, fixed, ext);
  std::copy(h1.rawptr(), h1.rawptr() + h1.size(), res.rawptr());

  return res;
}

} // namespace neml
