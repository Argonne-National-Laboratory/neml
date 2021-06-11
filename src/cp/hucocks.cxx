#include "hucocks.h"

namespace neml {

HuCocksPrecipitationModel::HuCocksPrecipitationModel(
    std::vector<std::shared_ptr<Interpolate>> c0, 
    std::vector<std::shared_ptr<Interpolate>> cp,
    std::vector<std::shared_ptr<Interpolate>> ceq,
    double am, double N0, double Vm, double chi, double D0, double Q0, 
    std::shared_ptr<Interpolate> Cf, double kboltz, double R,
    double Na, size_t rate, double f_init, double r_init, double N_init) :
      c0_(c0), cp_(cp), ceq_(ceq), am_(am), N0_(N0), Vm_(Vm), chi_(chi), 
      D0_(D0), Q0_(Q0), Cf_(Cf), kboltz_(kboltz), R_(R), Na_(Na),
      rate_(rate), f_init_(f_init), r_init_(r_init), N_init_(N_init),
      vm_(Vm / Na), varnames_({"f", "r", "N"})
{
}

std::string HuCocksPrecipitationModel::type()
{
  return "HuCocksPrecipitationModel";
}

std::unique_ptr<NEMLObject> HuCocksPrecipitationModel::initialize(
    ParameterSet & params)
{
  return neml::make_unique<HuCocksPrecipitationModel>(
      params.get_object_parameter_vector<Interpolate>("c0"),
      params.get_object_parameter_vector<Interpolate>("cp"),
      params.get_object_parameter_vector<Interpolate>("ceq"),
      params.get_parameter<double>("am"),
      params.get_parameter<double>("N0"),
      params.get_parameter<double>("Vm"),
      params.get_parameter<double>("chi"),
      params.get_parameter<double>("D0"),
      params.get_parameter<double>("Q0"),
      params.get_object_parameter<Interpolate>("Cf"),
      params.get_parameter<double>("kboltz"),
      params.get_parameter<double>("R"),
      params.get_parameter<double>("Na"),
      params.get_parameter<size_t>("rate"),
      params.get_parameter<double>("f_init"),
      params.get_parameter<double>("r_init"),
      params.get_parameter<double>("N_init"));
}

ParameterSet HuCocksPrecipitationModel::parameters()
{
  ParameterSet pset(HuCocksPrecipitationModel::type());

  pset.add_parameter<std::vector<NEMLObject>>("c0");
  pset.add_parameter<std::vector<NEMLObject>>("cp");
  pset.add_parameter<std::vector<NEMLObject>>("ceq");
  pset.add_parameter<double>("am");
  pset.add_parameter<double>("N0");
  pset.add_parameter<double>("Vm");
  pset.add_parameter<double>("chi");
  pset.add_parameter<double>("D0");
  pset.add_parameter<double>("Q0");
  pset.add_parameter<NEMLObject>("Cf");
  
  pset.add_optional_parameter<double>("kboltz", 1.3806485e-23);
  pset.add_optional_parameter<double>("R", 8.31462);
  pset.add_optional_parameter<double>("Na", 6.02e23);
  pset.add_optional_parameter<size_t>("rate", size_t(0));
  pset.add_optional_parameter<double>("f_init", 0);
  pset.add_optional_parameter<double>("r_init", 1e-20);
  pset.add_optional_parameter<double>("N_init", 1.0e-6);

  return pset;
}

std::vector<std::string> HuCocksPrecipitationModel::varnames() const
{
  return varnames_;
}

void HuCocksPrecipitationModel::set_varnames(std::vector<std::string> vars)
{
  if (vars.size() != varnames_.size())
    throw std::logic_error("New and old varname sizes do not match");
  varnames_ = vars;
}

void HuCocksPrecipitationModel::populate_history(History & history) const
{
  for (auto vn : varnames_)
    history.add<double>(vn);
}

void HuCocksPrecipitationModel::init_history(History & history) const
{
  history.get<double>(varnames_[0]) = f_init_;
  history.get<double>(varnames_[1]) = r_init_;
  history.get<double>(varnames_[2]) = N_init_;
}

double HuCocksPrecipitationModel::f_rate(const History & history, double T) const
{
  // double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);

  return 4.0 * M_PI / 3.0 * (N_rate(history, T) * std::pow(r, 3.0) + N * 3.0 *
                           std::pow(r, 2.0) * r_rate(history, T));
}

double HuCocksPrecipitationModel::df_df(const History & history, double T) const
{
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);

  return 4.0 * M_PI / 3.0 * (dN_df(history, T) * std::pow(r, 3.0) + N * 3.0 *
                           std::pow(r, 2.0) * dr_df(history, T));
}

double HuCocksPrecipitationModel::df_dr(const History & history, double T) const
{
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);

  return 4.0 * M_PI / 3.0 * (dN_dr(history, T) * std::pow(r, 3.0) +
                             3.0*N_rate(history, T) * std::pow(r, 2.0) + 6.0* N
                             * r * r_rate(history, T) + N * 3.0 * std::pow(r,
                                                                           2.0)
                             * dr_dr(history, T));
}

double HuCocksPrecipitationModel::df_dN(const History & history, double T) const
{
  // double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);

  return 4.0 * M_PI / 3.0 * (dN_dN(history, T) * std::pow(r, 3.0) + 3.0 *
                             std::pow(r, 2.0) * r_rate(history, T) + N * 3.0 *
                             std::pow(r, 2.0) * dr_dN(history, T));
}

double HuCocksPrecipitationModel::r_rate(const History & history, double T) const
{
  double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);

  auto ci = c(f, T);
  double D = D_(T);

  if (nucleation_(ci, T)) {
    double Gvi = Gv(f, T);
    double rc = -2.0 * chi_ / Gvi;
    return D / r * (ci[rate_] - ceq_[rate_]->value(T)) / (cp_[rate_]->value(T) -
                                                         ceq_[rate_]->value(T))
        + N_rate(history, T) / N * (rc - r);
  }
  else {
    double K = Cf_->value(T) * 8.0 * chi_ * Vm_ * D * ci[rate_] / (9.0 * R_ * T);
    return K / (3.0 * std::pow(r, 2.0));
  }
}

double HuCocksPrecipitationModel::dr_df(const History & history, double T) const
{
  double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);

  auto ci = c(f, T);
  auto dci = dc_df(f, T);
  double D = D_(T);

  if (nucleation_(ci, T)) {
    double Gvi = Gv(f, T);
    double dGvi = dG_df(f, T);
    double rc = -2.0 * chi_ / Gvi;
    double drc = 2.0 * chi_ / (Gvi * Gvi) * dGvi;
    return D / r * dci[rate_] / (cp_[rate_]->value(T) - ceq_[rate_]->value(T))
        + dN_df(history, T) / N * (rc - r) + N_rate(history, T) / N * drc;
  }
  else {
    double dK = Cf_->value(T) * 8.0 * chi_ * Vm_ * D * dci[rate_] / (9.0 * R_ * T);
    return dK / (3.0 * std::pow(r, 2.0));
  }
}

double HuCocksPrecipitationModel::dr_dr(const History & history, double T) const
{
  double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);

  auto ci = c(f, T);
  double D = D_(T);

  if (nucleation_(ci, T)) {
    double Gvi = Gv(f, T);
    double rc = -2.0 * chi_ / Gvi;
    return -D / (r*r) * (ci[rate_] - ceq_[rate_]->value(T)) / (cp_[rate_]->value(T) -
                                                         ceq_[rate_]->value(T))
        + dN_dr(history, T) / N * (rc - r) - N_rate(history,T) / N;
  }
  else {
    double K = Cf_->value(T) * 8.0 * chi_ * Vm_ * D * ci[rate_] / (9.0 * R_ * T);
    return -2 * K / (3.0 * std::pow(r, 3.0));
  }
}

double HuCocksPrecipitationModel::dr_dN(const History & history, double T) const
{
  double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);

  auto ci = c(f, T);

  if (nucleation_(ci, T)) {
    double Gvi = Gv(f, T);
    double rc = -2.0 * chi_ / Gvi;
    return dN_dN(history, T) / N * (rc - r) - N_rate(history,T)/(N*N) * (rc - r);
  }
  else {
    return 0;
  }
}

double HuCocksPrecipitationModel::N_rate(const History & history, double T) const
{
  double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);
  
  auto ci = c(f, T);
  double D = D_(T);

  if (nucleation_(ci, T)) {
    double Gvi = Gv(f, T);
    double Gstar = 16 * M_PI * std::pow(chi_, 3.0) / (3.0 * std::pow(Gvi, 2.0));
    double ZB = 2.0 * vm_ * D * ci[rate_] / std::pow(am_, 4.0) * std::sqrt(chi_ /
                                                                          (kboltz_
                                                                           * T));
    return N0_ * ZB * std::exp(-Gstar / (kboltz_ * T));
  }
  else {
    return -3.0 * N / r * r_rate(history, T);
  }
}

double HuCocksPrecipitationModel::dN_df(const History & history, double T) const
{
  double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);
  
  auto ci = c(f, T);
  auto dci = dc_df(f, T);
  double D = D_(T);

  if (nucleation_(ci, T)) {
    double Gvi = Gv(f, T);
    double dGvi = dG_df(f, T);
    double Gstar = 16 * M_PI * std::pow(chi_, 3.0) / (3.0 * std::pow(Gvi, 2.0));
    double dGstar = -32 * M_PI * std::pow(chi_, 3.0) / (3.0 * std::pow(Gvi,
                                                                       3.0)) *
        dGvi;

    double ZB = 2.0 * vm_ * D * ci[rate_] / std::pow(am_, 4.0) * std::sqrt(chi_ /
                                                                          (kboltz_
                                                                           * T));
    double dZB = 2.0 * vm_ * D * dci[rate_] / std::pow(am_, 4.0) * std::sqrt(chi_ /
                                                                          (kboltz_
                                                                           * T));

    return N0_ * dZB * std::exp(-Gstar / (kboltz_ * T)) - N0_ * ZB *
        std::exp(-Gstar/ (kboltz_*T)) * dGstar / (kboltz_ * T);
  }
  else {
    return -3.0 * N / r * dr_df(history, T);
  }
}

double HuCocksPrecipitationModel::dN_dr(const History & history, double T) const
{
  double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);
  
  auto ci = c(f, T);

  if (nucleation_(ci, T)) {
    return 0;
  }
  else {
    return -3.0 * N / r * dr_dr(history, T) + 3.0 * N / (r * r) * r_rate(history, T);
  }
}

double HuCocksPrecipitationModel::dN_dN(const History & history, double T) const
{
  double f = history.get<double>(varnames_[0]);
  double r = history.get<double>(varnames_[1]);
  double N = history.get<double>(varnames_[2]);
  
  auto ci = c(f, T);

  if (nucleation_(ci, T)) {
    return 0;
  }
  else {
    return -3.0 * N / r * dr_dN(history, T) - 3.0 / r * r_rate(history, T);
  }
}

size_t HuCocksPrecipitationModel::nspecies() const
{
  return c0_.size();
}

std::vector<double> HuCocksPrecipitationModel::c(double f, double T) const
{
  std::vector<double> cs(nspecies());
  for (size_t i = 0; i < nspecies(); i++)
    cs[i] = (c0_[i]->value(T) - f * cp_[i]->value(T)) / (1.0 - f);
  return cs;
}

std::vector<double> HuCocksPrecipitationModel::dc_df(double f, double T) const
{
  std::vector<double> cs(nspecies());
  for (size_t i = 0; i < nspecies(); i++)
    cs[i] = (c0_[i]->value(T) - cp_[i]->value(T)) / std::pow(1.0 - f, 2.0);
  return cs;
}

double HuCocksPrecipitationModel::D_(double T) const
{
  return D0_ * std::exp(-Q0_ / (R_ * T));
}

double HuCocksPrecipitationModel::Gv(double f, double T) const
{
  // Get the current concentrations
  auto ci = c(f, T);

  // Get the effective concentration
  double c_eff = 1;
  for (auto cii : ci)
    c_eff *= cii;
  
  // Get the effective equilibrium concentration
  double ceq_eff = 1;
  for (auto ci : ceq_)
    ceq_eff *= ci->value(T);

  return -kboltz_ * T / vm_ * std::log(c_eff / ceq_eff);
}

double HuCocksPrecipitationModel::dG_df(double f, double T) const
{
  // Get the current concentrations
  auto ci = c(f, T);

  // Get the derivative of the current concentrations
  auto dci = dc_df(f, T);

  // Get the effective concentration
  double c_eff = 1;
  for (auto cii : ci)
    c_eff *= cii;
  
  // Get the derivative of the effective concentration
  double dc_eff = 0.0;
  for (size_t i = 0; i < nspecies(); i++) {
    double a = 1.0;
    for (size_t j = 0; j < nspecies(); j++) {
      if (i == j)
        a *= dci[j];
      else
        a *= ci[j];
    }
    dc_eff += a;
  }
  
  // Get the effective equilibrium concentration
  double ceq_eff = 1;
  for (auto ci : ceq_)
    ceq_eff *= ci->value(T);

  return -kboltz_ * T / vm_ * dc_eff / c_eff;
}

bool HuCocksPrecipitationModel::nucleation_(const std::vector<double> & c,
                                            double T) const
{
  for (size_t i = 0; i < nspecies(); i++) {
    if (c[i] <= ceq_[i]->value(T)) return false;
  }
  return true;
}

DislocationSpacingHardening::DislocationSpacingHardening(std::shared_ptr<Interpolate> J1, 
                                                         std::shared_ptr<Interpolate> J2,
                                                         std::shared_ptr<Interpolate> K, 
                                                         double L0,
                                                         double a, double b,
                                                         std::shared_ptr<Interpolate> G,
                                                         std::shared_ptr<Lattice> L,
                                                         std::string varprefix) :
    J1_(J1), J2_(J2), K_(K), L0_(L0), a_(a), b_(b), G_(G), L_(L),
    varprefix_(varprefix)
{
  varnames_.resize(L_->ntotal());
  for (size_t i = 0; i < size(); i++)
    varnames_[i] = varprefix_+"_"+std::to_string(i);

  init_cache_();
}

std::string DislocationSpacingHardening::type()
{
  return "DislocationSpacingHardening";
}

std::unique_ptr<NEMLObject> DislocationSpacingHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<DislocationSpacingHardening>(
      params.get_object_parameter<Interpolate>("J1"),
      params.get_object_parameter<Interpolate>("J2"),
      params.get_object_parameter<Interpolate>("K"),
      params.get_parameter<double>("L0"),
      params.get_parameter<double>("a"),
      params.get_parameter<double>("b"),
      params.get_object_parameter<Interpolate>("G"),
      params.get_object_parameter<Lattice>("L"),
      params.get_parameter<std::string>("varprefix"));
}

ParameterSet DislocationSpacingHardening::parameters()
{
  ParameterSet pset(DislocationSpacingHardening::type());

  pset.add_parameter<NEMLObject>("J1");
  pset.add_parameter<NEMLObject>("J2");
  pset.add_parameter<NEMLObject>("K");
  pset.add_parameter<double>("L0");
  pset.add_parameter<double>("a");
  pset.add_parameter<double>("b");
  pset.add_parameter<NEMLObject>("G");
  pset.add_parameter<NEMLObject>("L");

  pset.add_optional_parameter<std::string>("varprefix", 
                                           std::string("spacing"));

  return pset;
}

std::vector<std::string> DislocationSpacingHardening::varnames() const
{
  return varnames_;
}

void DislocationSpacingHardening::set_varnames(std::vector<std::string> vars)
{
  varnames_ = vars;
  init_cache_();
}

void DislocationSpacingHardening::populate_history(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void DislocationSpacingHardening::init_history(History & history) const
{
  for (auto vn : varnames_) {
    history.get<double>(vn) = L0_;
  }
}

double DislocationSpacingHardening::hist_to_tau(size_t g, size_t i, 
                                           const History & history,
                                           Lattice & L,
                                           double T, const History & fixed) const
{
  return a_ * G_->value(T) * b_ / history.get<double>(varnames_[L.flat(g,i)]);
}

History DislocationSpacingHardening::d_hist_to_tau(size_t g, size_t i, 
                                              const History & history,
                                              Lattice & L,
                                              double T, 
                                              const History & fixed) const
{
  History res = cache(CacheType::DOUBLE);
  // This works because the above zeros out the vector
  res.get<double>(varnames_[L.flat(g,i)]) = -a_ * G_->value(T) * b_ / 
      std::pow(history.get<double>(varnames_[L.flat(g,i)]), 2.0);
  return res;
}

History DislocationSpacingHardening::hist(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history, 
                                     Lattice & L, double T, const SlipRule & R, 
                                     const History & fixed) const
{

  // Vector of results
  History res = blank_hist().zero();
  
  for (size_t i = 0; i < size(); i++) {
    double Li = history.get<double>(varnames_[i]);
    double Li3 = std::pow(Li, 3);
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t k = 0; k < L.nslip(g); k++) {
        size_t j = L.flat(g, k);
        double c;
        if (i == j)
          c = J1_->value(T);
        else
          c = J2_->value(T);
        res.get<double>(varnames_[i]) -= Li3 * c * 
            std::fabs(R.slip(g, k, stress, Q, history, L, T, fixed));
      }
    }
    res.get<double>(varnames_[i]) += K_->value(T) / Li3;
  }


  return res;
}

History DislocationSpacingHardening::d_hist_d_s(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history,
                                           Lattice & L, double T, 
                                           const SlipRule & R,
                                           const History & fixed) const
{
  History res = blank_hist().derivative<Symmetric>().zero();

  for (size_t i = 0; i < size(); i++) {
    double Li = history.get<double>(varnames_[i]);
    double Li3 = std::pow(Li, 3);
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t k = 0; k < L.nslip(g); k++) {
        size_t j = L.flat(g, k);
        double c;
        if (i == j)
          c = J1_->value(T);
        else
          c = J2_->value(T);
        
        double si = R.slip(g, k, stress, Q, history, L, T, fixed);
        double sgn = std::copysign(1.0, si);
        Symmetric dsi = R.d_slip_d_s(g, k, stress, Q, history, L, T, fixed);
        
        res.get<Symmetric>(varnames_[i]) -= Li3 * c * sgn * dsi;
      }
    }
  }

  return res;
}

History DislocationSpacingHardening::d_hist_d_h(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history, 
                                           Lattice & L,
                                           double T, const SlipRule & R, 
                                           const History & fixed) const
{
  auto res = blank_hist().derivative<History>().zero();

  // There is both an effect on the diagonal caused by the Ldi in the equation
  // and an effect from the slip rule
  for (size_t i = 0; i < size(); i++) {
    double Li = history.get<double>(varnames_[i]);
    double Li3 = std::pow(Li, 3);
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t k = 0; k < L.nslip(g); k++) {
        size_t j = L.flat(g, k);
        double c;
        if (i == j)
          c = J1_->value(T);
        else
          c = J2_->value(T);
        
        double si = R.slip(g, k, stress, Q, history, L, T, fixed);
        double sgn = std::copysign(1.0, si);
        History dsi = R.d_slip_d_h(g, k, stress, Q, history, L, T, fixed);

        res.get<double>(varnames_[i]+"_"+varnames_[i]) -= 
            3.0 * std::pow(Li,2.0) * c * std::fabs(si);
        
        for (auto vn : varnames_)
          res.get<double>(varnames_[i]+"_"+vn) -=
              Li3 * c * sgn * dsi.get<double>(vn);
      }
    }
    res.get<double>(varnames_[i]+"_"+varnames_[i]) -= 3.0*K_->value(T) /
        std::pow(Li,4.0);    
  }

  return res;
}

History DislocationSpacingHardening::d_hist_d_h_ext(const Symmetric & stress, 
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & L, double T, const SlipRule & R,
                                               const History & fixed, 
                                               std::vector<std::string> ext) const
{
  History res = blank_hist().history_derivative(history.subset(ext)).zero();

  for (size_t i = 0; i < size(); i++) {
    double Li = history.get<double>(varnames_[i]);
    double Li3 = std::pow(Li, 3);
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t k = 0; k < L.nslip(g); k++) {
        size_t j = L.flat(g, k);
        double c;
        if (i == j)
          c = J1_->value(T);
        else
          c = J2_->value(T);
        
        double si = R.slip(g, k, stress, Q, history, L, T, fixed);
        double sgn = std::copysign(1.0, si);
        History dsi = R.d_slip_d_h(g, k, stress, Q, history, L, T, fixed);

        for (auto vn : ext)
          if (dsi.contains(vn)) {
            res.get<double>(varnames_[i]+"_"+vn) -=
                Li3 * c * sgn * dsi.get<double>(vn);
          }
      }
    }
  }

  return res;
}

size_t DislocationSpacingHardening::size() const
{
  return varnames_.size();
}

}
