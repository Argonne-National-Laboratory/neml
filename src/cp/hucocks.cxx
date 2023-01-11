#include "cp/hucocks.h"

namespace neml {

HuCocksPrecipitationModel::HuCocksPrecipitationModel(
    ParameterSet & params) :
      HistoryNEMLObject(params),
      c0_(params.get_object_parameter_vector<Interpolate>("c0")), 
      cp_(params.get_object_parameter_vector<Interpolate>("cp")), 
      ceq_(params.get_object_parameter_vector<Interpolate>("ceq")),
      am_(params.get_parameter<double>("am")), 
      N0_(params.get_parameter<double>("N0")), 
      Vm_(params.get_parameter<double>("Vm")), 
      chi_(params.get_parameter<double>("chi")), 
      D0_(params.get_parameter<double>("D0")), 
      Q0_(params.get_parameter<double>("Q0")), 
      Cf_(params.get_object_parameter<Interpolate>("Cf")), 
      kboltz_(params.get_parameter<double>("kboltz")), 
      R_(params.get_parameter<double>("R")), 
      Na_(params.get_parameter<double>("Na")),
      rate_(params.get_parameter<size_t>("rate")), 
      f_init_(params.get_parameter<double>("f_init")), 
      r_init_(params.get_parameter<double>("r_init")), 
      N_init_(params.get_parameter<double>("N_init")),
      fs_(params.get_parameter<double>("fs")), 
      rs_(params.get_parameter<double>("rs")), 
      Ns_(params.get_parameter<double>("Ns")), 
      w_(params.get_parameter<double>("w")),
      vm_(Vm_ / Na_), 
      varnames_({"f", "r", "N"})
{

}

std::string HuCocksPrecipitationModel::type()
{
  return "HuCocksPrecipitationModel";
}

std::unique_ptr<NEMLObject> HuCocksPrecipitationModel::initialize(
    ParameterSet & params)
{
  return neml::make_unique<HuCocksPrecipitationModel>(params);
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
  pset.add_optional_parameter<double>("f_init", 4.18879e-16);
  pset.add_optional_parameter<double>("r_init", 1.0e-9);
  pset.add_optional_parameter<double>("N_init", 1.0e11);
  pset.add_optional_parameter<double>("fs", 0.1);
  pset.add_optional_parameter<double>("rs", 1.0e-9);
  pset.add_optional_parameter<double>("Ns", 1.0e12);
  pset.add_optional_parameter<double>("w", 1.0);

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

void HuCocksPrecipitationModel::populate_hist(History & history) const
{
  for (auto vn : varnames_)
    history.add<double>(vn);
}

void HuCocksPrecipitationModel::init_hist(History & history) const
{
  history.get<double>(varnames_[0]) = f_init_ / fs_;
  history.get<double>(varnames_[1]) = r_init_ / rs_;
  history.get<double>(varnames_[2]) = N_init_ / Ns_;
}

std::vector<double> HuCocksPrecipitationModel::rate(const History & history,
                                                    double T) const
{
  double fi = f(history);
  double ri = r(history);
  double Ni = N(history);

  return {f_rate(fi, ri, Ni, T) / fs_, r_rate(fi, ri, Ni, T) / rs_, 
    N_rate(fi, ri, Ni, T) / Ns_};
}

std::vector<std::vector<double>> HuCocksPrecipitationModel::jac(
    const History & history,double T) const
{
  double fi = f(history);
  double ri = r(history);
  double Ni = N(history);

  return
  {
    {df_df(fi, ri, Ni, T) / fs_ * fs_, df_dr(fi, ri, Ni, T) / fs_ * rs_, df_dN(fi, ri, Ni, T) / fs_ * Ns_},
    {dr_df(fi, ri, Ni, T) / rs_ * fs_, dr_dr(fi, ri, Ni, T) / rs_ * rs_, dr_dN(fi, ri, Ni, T) / rs_ * Ns_},
    {dN_df(fi, ri, Ni, T) / Ns_ * fs_, dN_dr(fi, ri, Ni, T) / Ns_ * rs_, dN_dN(fi, ri, Ni, T )/ Ns_ * Ns_}
  };
}

double HuCocksPrecipitationModel::f(const History & history) const
{
  return history.get<double>(varnames_[0]) * fs_;
}

double HuCocksPrecipitationModel::r(const History & history) const
{
  return history.get<double>(varnames_[1]) * rs_;
}

double HuCocksPrecipitationModel::N(const History & history) const
{
  return history.get<double>(varnames_[2]) * Ns_;
}

double HuCocksPrecipitationModel::f_rate(double f, double r, double N, double T) const
{
  return  4.0 * M_PI / 3.0 * (N_rate(f, r, N, T) * std::pow(r, 3.0) + N * 3.0 *
                           std::pow(r, 2.0) * r_rate(f, r, N, T));
}

double HuCocksPrecipitationModel::df_df(double f, double r, double N, double T) const
{
  return 4.0 * M_PI / 3.0 * (dN_df(f, r, N, T) * std::pow(r, 3.0) + N * 3.0 *
                           std::pow(r, 2.0) * dr_df(f, r, N, T));
}

double HuCocksPrecipitationModel::df_dr(double f, double r, double N, double T) const
{
  return 4.0 * M_PI / 3.0 * (dN_dr(f, r, N, T) * std::pow(r, 3.0) +
                             3.0*N_rate(f, r, N, T) * std::pow(r, 2.0) + 6.0* N
                             * r * r_rate(f, r, N, T) + N * 3.0 * std::pow(r,
                                                                           2.0)
                             * dr_dr(f, r, N, T));
}

double HuCocksPrecipitationModel::df_dN(double f, double r, double N, double T) const
{
  return 4.0 * M_PI / 3.0 * (dN_dN(f, r, N, T) * std::pow(r, 3.0) + 3.0 *
                             std::pow(r, 2.0) * r_rate(f, r, N, T) + N * 3.0 *
                             std::pow(r, 2.0) * dr_dN(f, r, N, T));
}

double HuCocksPrecipitationModel::r_rate(double f, double r, double N, double T) const
{
  double s, ds;
  sfn_(f, T, s, ds);

  return (1.0-s) * r_rate_nucleation_(f, r, N, T) + s * r_rate_ripening_(f, r, N,
                                                                       T);
}

double HuCocksPrecipitationModel::dr_df(double f, double r, double N, double T) const
{
  double s, ds;
  sfn_(f, T, s, ds);

  return (1.0-s) * dr_df_nucleation_(f, r, N, T) + s * dr_df_ripening_(f, r, N, T) -
      ds * r_rate_nucleation_(f, r, N, T) + ds * r_rate_ripening_(f, r, N, T);
}

double HuCocksPrecipitationModel::dr_dr(double f, double r, double N, double T) const
{
  double s, ds;
  sfn_(f, T, s, ds);

  return (1.0-s) * dr_dr_nucleation_(f, r, N, T) + s * dr_dr_ripening_(f, r, N, T);
}

double HuCocksPrecipitationModel::dr_dN(double f, double r, double N, double T) const
{
  double s, ds;
  sfn_(f, T, s, ds);

  return (1.0-s) * dr_dN_nucleation_(f, r, N, T) + s * dr_dN_ripening_(f, r, N,
                                                                     T);
}

double HuCocksPrecipitationModel::r_rate_nucleation_(double f, double r, double N, double T) const
{
  auto ci = c(f, T);
  double D = D_(T);

  double Gvi = Gv(f, T);
  double rc = -2.0 * chi_ / Gvi;
  return D / r * (ci[rate_] - ceq_[rate_]->value(T)) / (cp_[rate_]->value(T) -
                                                       ceq_[rate_]->value(T))
      + N_rate_nucleation_(f, r, N, T) / N * (rc - r);

}

double HuCocksPrecipitationModel::dr_df_nucleation_(double f, double r, double N, double T) const
{
  auto ci = c(f, T);
  auto dci = dc_df(f, T);
  double D = D_(T);

  double Gvi = Gv(f, T);
  double dGvi = dG_df(f, T);
  double rc = -2.0 * chi_ / Gvi;
  double drc = 2.0 * chi_ / (Gvi * Gvi) * dGvi;
  return D / r * dci[rate_] / (cp_[rate_]->value(T) - ceq_[rate_]->value(T))
      + dN_df_nucleation_(f, r, N, T) / N * (rc - r) + N_rate_nucleation_(f, r, N, T) / N * drc;

}

double HuCocksPrecipitationModel::dr_dr_nucleation_(double f, double r, double N, double T) const
{
  auto ci = c(f, T);
  double D = D_(T);

  double Gvi = Gv(f, T);
  double rc = -2.0 * chi_ / Gvi;
  return -D / (r*r) * (ci[rate_] - ceq_[rate_]->value(T)) / (cp_[rate_]->value(T) -
                                                       ceq_[rate_]->value(T))
      + dN_dr_nucleation_(f, r, N, T) / N * (rc - r) - N_rate_nucleation_(f, r, N,T) / N;
}

double HuCocksPrecipitationModel::dr_dN_nucleation_(double f, double r, double N, double T) const
{
  auto ci = c(f, T);

  double Gvi = Gv(f, T);
  double rc = -2.0 * chi_ / Gvi;
  return dN_dN_nucleation_(f, r, N, T) / N * (rc - r) - N_rate_nucleation_(f, r, N,T)/(N*N) * (rc - r);
}

double HuCocksPrecipitationModel::r_rate_ripening_(double f, double r, double N, double T) const
{
  auto ci = c(f, T);
  double D = D_(T);

  double K = Cf_->value(T) * 8.0 * chi_ * Vm_ * D * ci[rate_] / (9.0 * R_ * T);
  return K / (3.0 * std::pow(r, 2.0));
}

double HuCocksPrecipitationModel::dr_df_ripening_(double f, double r, double N, double T) const
{
  auto dci = dc_df(f, T);
  double D = D_(T);

  double dK = Cf_->value(T) * 8.0 * chi_ * Vm_ * D * dci[rate_] / (9.0 * R_ * T);
  return dK / (3.0 * std::pow(r, 2.0));
}

double HuCocksPrecipitationModel::dr_dr_ripening_(double f, double r, double N, double T) const
{
  auto ci = c(f, T);
  double D = D_(T);

  double K = Cf_->value(T) * 8.0 * chi_ * Vm_ * D * ci[rate_] / (9.0 * R_ * T);
  return -2 * K / (3.0 * std::pow(r, 3.0));
}

double HuCocksPrecipitationModel::dr_dN_ripening_(double f, double r, double N, double T) const
{
  return 0;
}

double HuCocksPrecipitationModel::N_rate(double f, double r, double N, double T) const
{
  double s, ds;
  sfn_(f, T, s, ds);

  return (1.0-s) * N_rate_nucleation_(f, r, N, T) + s * N_rate_ripening_(f, r, N,
                                                                       T);
}

double HuCocksPrecipitationModel::dN_df(double f, double r, double N, double T) const
{
  double s, ds;
  sfn_(f, T, s, ds);

  return (1.0-s) * dN_df_nucleation_(f, r, N, T) + s * dN_df_ripening_(f, r, N, T) -
      ds * N_rate_nucleation_(f, r, N, T) + ds * N_rate_ripening_(f, r, N, T);
}

double HuCocksPrecipitationModel::dN_dr(double f, double r, double N, double T) const
{
  double s, ds;
  sfn_(f, T, s, ds);

  return (1.0-s) * dN_dr_nucleation_(f, r, N, T) + s * dN_dr_ripening_(f, r, N, T);
}

double HuCocksPrecipitationModel::dN_dN(double f, double r, double N, double T) const
{
  double s, ds;
  sfn_(f, T, s, ds);

  return (1.0-s) * dN_dN_nucleation_(f, r, N, T) + s * dN_dN_ripening_(f, r, N,
                                                                     T);
}

double HuCocksPrecipitationModel::N_rate_nucleation_(double f, double r, double N, double T) const
{
  auto ci = c(f, T);
  double D = D_(T);

  double Gvi = Gv(f, T);
  double Gstar = 16 * M_PI * std::pow(chi_, 3.0) / (3.0 * std::pow(Gvi, 2.0))
      * w_;
  double ZB = 2.0 * vm_ * D * ci[rate_] / std::pow(am_, 4.0) * std::sqrt(chi_ /
                                                                        (kboltz_
                                                                         * T));
  return N0_ * ZB * std::exp(-Gstar / (kboltz_ * T));

}

double HuCocksPrecipitationModel::dN_df_nucleation_(double f, double r, double N, double T) const
{
  auto ci = c(f, T);
  auto dci = dc_df(f, T);
  double D = D_(T);

  double Gvi = Gv(f, T);
  double dGvi = dG_df(f, T);
  double Gstar = 16 * M_PI * std::pow(chi_, 3.0) / (3.0 * std::pow(Gvi, 2.0))
      * w_;
  double dGstar = -32 * M_PI * std::pow(chi_, 3.0) / (3.0 * std::pow(Gvi,
                                                                     3.0)) *
      dGvi * w_;

  double ZB = 2.0 * vm_ * D * ci[rate_] / std::pow(am_, 4.0) * std::sqrt(chi_ /
                                                                        (kboltz_
                                                                         * T));
  double dZB = 2.0 * vm_ * D * dci[rate_] / std::pow(am_, 4.0) * std::sqrt(chi_ /
                                                                        (kboltz_
                                                                         * T));

  return N0_ * dZB * std::exp(-Gstar / (kboltz_ * T)) - N0_ * ZB *
      std::exp(-Gstar/ (kboltz_*T)) * dGstar / (kboltz_ * T);
}

double HuCocksPrecipitationModel::dN_dr_nucleation_(double f, double r, double N, double T) const
{
  return 0;
}

double HuCocksPrecipitationModel::dN_dN_nucleation_(double f, double r, double N, double T) const
{
  return 0;
}

double HuCocksPrecipitationModel::N_rate_ripening_(double f, double r, double N, double T) const
{
  return -3.0 * N / r * r_rate_ripening_(f, r, N, T);
}

double HuCocksPrecipitationModel::dN_df_ripening_(double f, double r, double N, double T) const
{
  return -3.0 * N / r * dr_df_ripening_(f, r, N, T);
}

double HuCocksPrecipitationModel::dN_dr_ripening_(double f, double r, double N, double T) const
{
  return -3.0 * N / r * dr_dr_ripening_(f, r, N, T) + 3.0 * N / (r * r) * r_rate_ripening_(f, r, N, T);
}

double HuCocksPrecipitationModel::dN_dN_ripening_(double f, double r, double N, double T) const
{
  return -3.0 * N / r * dr_dN_ripening_(f, r, N, T) - 3.0 / r * r_rate_ripening_(f, r, N, T);
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
  
  return -kboltz_ * T / vm_ * dc_eff / c_eff;
}

double HuCocksPrecipitationModel::vm() const
{
  return vm_;
}

void HuCocksPrecipitationModel::sfn_(double f, double T, double & val, double
                                       & dval) const
{
  auto cs = c(f, T);
  auto dcs = dc_df(f, T);
  val = 0.0;
  for (size_t i = 0; i < nspecies(); i++) {
    double xi = (cs[i] - c0_[i]->value(T)) / (ceq_[i]->value(T) -
                                              c0_[i]->value(T));
    if (xi > val) {
      val = xi;
      dval = dcs[i] / (ceq_[i]->value(T) - c0_[i]->value(T));
    }
  }
  
  if (val < 0.0) {
    val = 0.0;
    dval = 0.0;
  }
  
  if (val > 1.0) {
    val = 1.0;
    dval = 0.0;
  }
}

bool HuCocksPrecipitationModel::nucleation_(const std::vector<double> & c,
                                            double T) const
{
  for (size_t i = 0; i < nspecies(); i++) {
    if (c[i] <= ceq_[i]->value(T)) return false;
  }
  return true;
}

DislocationSpacingHardening::DislocationSpacingHardening(ParameterSet & params) :
    SlipHardening(params),
    J1_(params.get_object_parameter<Interpolate>("J1")),
    J2_(params.get_object_parameter<Interpolate>("J2")),
    K_(params.get_object_parameter<Interpolate>("K")),
    L0_(params.get_parameter<double>("L0")),
    a_(params.get_parameter<double>("a")),
    b_(params.get_parameter<double>("b")),
    G_(params.get_object_parameter<Interpolate>("G")),
    L_(params.get_object_parameter<Lattice>("L")),
    varprefix_(params.get_parameter<std::string>("varprefix"))
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
  return neml::make_unique<DislocationSpacingHardening>(params);
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

void DislocationSpacingHardening::populate_hist(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void DislocationSpacingHardening::init_hist(History & history) const
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
  auto res = blank_hist().history_derivative(history).zero();

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
        
        for (auto vn : dsi.get_order())
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

HuCocksHardening::HuCocksHardening(ParameterSet & params) :
    SlipHardening(params),
    dmodel_(params.get_object_parameter<SlipHardening>("dmodel")),
    pmodels_(params.get_object_parameter_vector<HuCocksPrecipitationModel>("pmodels")),
    ap_(params.get_parameter<double>("ap")),
    ac_(params.get_parameter<double>("ac")),
    b_(params.get_parameter<double>("b")),
    G_(params.get_object_parameter<Interpolate>("G"))
{
  // Alter the varnames of the pmodels...
  size_t i = 0;
  for (auto pmodel : pmodels_) {
    std::vector<std::string> new_varnames;
    for (auto vn : pmodel->varnames()) {
      new_varnames.push_back(vn + "_" + std::to_string(i));
    }
    pnames_.push_back(new_varnames);
    pmodel->set_varnames(new_varnames);
    i++;
  }
  init_cache_();
}

std::string HuCocksHardening::type()
{
  return "HuCocksHardening";
}

std::unique_ptr<NEMLObject> HuCocksHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<HuCocksHardening>(params);
}

ParameterSet HuCocksHardening::parameters()
{
  ParameterSet pset(HuCocksHardening::type());

  pset.add_parameter<NEMLObject>("dmodel");
  pset.add_parameter<std::vector<NEMLObject>>("pmodels");
  pset.add_parameter<double>("ap");
  pset.add_parameter<double>("ac");
  pset.add_parameter<double>("b");
  pset.add_parameter<NEMLObject>("G");

  return pset;
}

std::vector<std::string> HuCocksHardening::varnames() const
{
  std::vector<std::string> names = dmodel_->varnames();
  for (auto pmodel : pmodels_) {
    auto pn = pmodel->varnames();
    names.insert(names.end(), pn.begin(), pn.end());
  }

  return names;
}

void HuCocksHardening::set_varnames(std::vector<std::string> vars)
{
  throw std::runtime_error("Cannot override varnames for HuCocksHardening");
}

void HuCocksHardening::populate_hist(History & history) const
{
  dmodel_->populate_hist(history);
  for (auto pmodel : pmodels_)
    pmodel->populate_hist(history);
}

void HuCocksHardening::init_hist(History & history) const
{
  history.zero();
  dmodel_->init_hist(history);
  for (auto pmodel : pmodels_)
    pmodel->init_hist(history);
}

double HuCocksHardening::hist_to_tau(size_t g, size_t i, 
                                           const History & history,
                                           Lattice & L,
                                           double T, const History & fixed) const
{
  double tau_d = dmodel_->hist_to_tau(g, i, history, L, T, fixed);
  double c = c_eff_(history, T);
  double NA = NA_eff_(history, T); 

  double tau_p = ap_ * G_->value(T) * b_ * std::sqrt(NA);
  double tau_c = ac_ * G_->value(T) * b_ * std::sqrt(c * b_);

  return std::sqrt(tau_d*tau_d + tau_p*tau_p) + tau_c;
}

History HuCocksHardening::d_hist_to_tau(size_t g, size_t i, 
                                              const History & history,
                                              Lattice & L,
                                              double T, 
                                              const History & fixed) const
{
  History res = cache(CacheType::DOUBLE).zero();

  // Commonly-used things
  double tau_d = dmodel_->hist_to_tau(g, i, history, L, T, fixed);
  double c = c_eff_(history, T);
  double NA = NA_eff_(history, T); 

  double tau_p = ap_ * G_->value(T) * b_ * std::sqrt(NA);

  // First block: dislocation terms
  History dd = dmodel_->d_hist_to_tau(g, i, history, L, T, fixed);
  dd.scalar_multiply(tau_d / std::sqrt(tau_p * tau_p + tau_d * tau_d));
  std::copy(dd.rawptr(), dd.rawptr()+dd.size(), res.rawptr());
  
  // For each precipitation model
  for (size_t i = 0; i < pmodels_.size(); i++) {
    // Second block: f
    auto dc = pmodels_[i]->dc_df(pmodels_[i]->f(history),T);
    for (auto dci: dc)
      res.get<double>(pnames_[i][0]) += (ac_ * G_->value(T) * b_ * b_ / (2.0 *
                                                                       std::sqrt(c
                                                                                 * b_))
          * dci / pmodels_[i]->vm())*pmodels_[i]->fs();

    // Third block: r
    res.get<double>(pnames_[i][1]) = (tau_p / std::sqrt(tau_p * tau_p + tau_d *
                                                       tau_d) * ap_ *
        G_->value(T) * b_ * pmodels_[i]->N(history) / std::sqrt(NA))
        * pmodels_[i]->rs();

    // Fourth block: N
    res.get<double>(pnames_[i][2]) = (tau_p / std::sqrt(tau_p * tau_p + tau_d *
                                                       tau_d) * ap_ *
        G_->value(T) * b_ * pmodels_[i]->r(history) / std::sqrt(NA)) *
        pmodels_[i]->Ns();
  }

  return res;
}

History HuCocksHardening::hist(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history, 
                                     Lattice & L, double T, const SlipRule & R, 
                                     const History & fixed) const
{

  // Vector of results
  History res = blank_hist().zero();
  
  // History is easy, we just concatenate the dislocation model hardening
  // followed by each precipitation hardening
  auto h1 = dmodel_->hist(stress, Q, history, L, T, R, fixed);
  std::copy(h1.rawptr(), h1.rawptr() + h1.size(), 
            res.start_loc(dmodel_->varnames()[0]));
  for (size_t i = 0; i < pmodels_.size(); i++) {
    auto rate = pmodels_[i]->rate(history, T);
    for (size_t j = 0; j < 3; j++)
      res.get<double>(pnames_[i][j]) = rate[j];
  }

  return res;
}

History HuCocksHardening::d_hist_d_s(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history,
                                           Lattice & L, double T, 
                                           const SlipRule & R,
                                           const History & fixed) const
{
  // This could be non-zero
  History res = dmodel_->d_hist_d_s(stress, Q, history, L, T, R, fixed);

  // These are zero
  for (size_t i = 0; i < pmodels_.size(); i++) {
    History pm;
    pmodels_[i]->populate_hist(pm);
    History zero = pm.derivative<Symmetric>().zero();
    res.add_union(zero);
  }

  return res;
}

History HuCocksHardening::d_hist_d_h(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history, 
                                           Lattice & L,
                                           double T, const SlipRule & R, 
                                           const History & fixed) const
{
  // Cache some useful objects
  std::vector<History> phists;
  for (auto pmodel : pmodels_) {
    History pm;
    pmodel->populate_hist(pm);
    phists.push_back(pm.zero());
  }

  // Start with the self-derivative
  auto res = dmodel_->d_hist_d_h(stress, Q, history, L, T, R, fixed);

  // Deal with both cross and actual terms
  for (size_t i = 0; i < pmodels_.size(); i++) {
    // dmodel/pmodel cross term
    //res.add_union(dmodel_->blank_hist().history_derivative(phists[i]).zero());
    // pmodel/dmodel cross term
    res.add_union(phists[i].history_derivative(dmodel_->blank_hist()).zero());
    // pmodel/other pmodel cross terms
    for (size_t j = 0; j < pmodels_.size(); j++) {
      if (i != j)
        res.add_union(phists[i].history_derivative(phists[j]).zero());
    }
    // actual non-zeros
    auto jac = pmodels_[i]->jac(history, T);
    for (size_t ii = 0; ii < 3; ii++) {
      for (size_t jj = 0; jj < 3; jj++) {
        res.add<double>(pnames_[i][ii]+"_"+pnames_[i][jj]);
        res.get<double>(pnames_[i][ii]+"_"+pnames_[i][jj]) = jac[ii][jj];
      }
    }
  }

  // Reorder...
  std::vector<std::string> order;
  for (auto n1 :  varnames()) {
    for (auto n2 : varnames()) {
      order.push_back(n1 + "_" + n2);
    }
  }
  res.reorder(order);

  return res;
}

History HuCocksHardening::d_hist_d_h_ext(const Symmetric & stress, 
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & L, double T, const SlipRule & R,
                                               const History & fixed, 
                                               std::vector<std::string> ext) const
{
  History res = blank_hist().history_derivative(history.subset(ext)).zero();

  // The only potential non-zero is from the dmodel
  History h1 = dmodel_->d_hist_d_h_ext(stress, Q, history, L, T, R, 
                                         fixed, ext);
  std::copy(h1.rawptr(), h1.rawptr() + h1.size(), res.rawptr());

  return res;
}

double HuCocksHardening::c_eff_(const History & history, double T) const
{
  double res = 0.0;
  for (size_t i = 0; i < pmodels_.size(); i++) {
    auto c = pmodels_[i]->c(pmodels_[i]->f(history), T);
    for (auto ci : c)
      res += ci / pmodels_[i]->vm();
  }
  return res;
}

double HuCocksHardening::NA_eff_(const History & history, double T) const
{
  double res = 0.0;
  for (size_t i = 0; i < pmodels_.size(); i++) {
    res += 2.0 * pmodels_[i]->r(history) * pmodels_[i]->N(history);
  }
  return res;
}

ArrheniusSlipRule::ArrheniusSlipRule(ParameterSet & params) :
    SlipStrengthSlipRule(params),
    g0_(params.get_parameter<double>("g0")),
    A_(params.get_parameter<double>("A")),
    B_(params.get_parameter<double>("B")),
    b_(params.get_parameter<double>("b")),
    a0_(params.get_parameter<double>("a0")),
    G0_(params.get_parameter<double>("G0")),
    k_(params.get_parameter<double>("k"))
{

}

std::string ArrheniusSlipRule::type()
{
  return "ArrheniusSlipRule";
}

std::unique_ptr<NEMLObject> ArrheniusSlipRule::initialize(
    ParameterSet & params)
{
  return neml::make_unique<ArrheniusSlipRule>(params);
}

ParameterSet ArrheniusSlipRule::parameters()
{
  ParameterSet pset(ArrheniusSlipRule::type());
  
  pset.add_parameter<NEMLObject>("resistance");
  pset.add_parameter<double>("g0");
  pset.add_parameter<double>("A");
  pset.add_parameter<double>("B");
  pset.add_parameter<double>("b");
  pset.add_parameter<double>("a0");
  pset.add_parameter<double>("G0");
  pset.add_optional_parameter<double>("k", 1.3806485e-23); 

  return pset;
}

double ArrheniusSlipRule::scalar_sslip(size_t g, size_t i, double tau, 
                                      double strength, double T) const
{
  if (tau == 0.0)
    return 0.0;

  double F0kT = a0_ * G0_ * std::pow(b_,3.0) / (k_ * T);

  return  g0_ * std::exp(-F0kT * std::pow(1.0 - std::pow(std::fabs(tau/strength),
                                                     A_), B_)) * 
      std::copysign(1.0, tau);
}

double ArrheniusSlipRule::scalar_d_sslip_dtau(size_t g, size_t i, double tau, 
                                             double strength, double T) const
{
  if (tau == 0.0)
    return 0.0;

  double F0kT = a0_ * G0_ * std::pow(b_,3.0) / (k_ * T);

  return std::fabs(-g0_ * A_*B_*F0kT*tau*std::pow(std::fabs(tau/strength),
                                  A_-2.0)*std::pow(1.0-std::pow(std::fabs(tau/strength),
                                                                A_), B_-1.0) *  
                                      std::exp(-F0kT * std::pow(1.0 -
                                                                std::pow(std::fabs(tau/strength),
                                                                         A_),
                                                                B_)) /
                                      std::pow(strength, 2.0));
}

double ArrheniusSlipRule::scalar_d_sslip_dstrength(size_t g, size_t i, 
                                                  double tau, 
                                                  double strength,
                                                  double T) const
{
  if (tau == 0.0)
    return 0.0;

  double F0kT = a0_ * G0_ * std::pow(b_,3.0) / (k_ * T);
  
  return  -g0_ * A_ * B_ * F0kT * std::pow(tau, 2.0) *
                       std::pow(std::fabs(tau/strength), A_-2.0) * std::pow(1.0
                                                                            - std::pow(std::fabs(tau/strength),
                                                                                       A_),
                                                                            B_-1.0)
                       * std::exp(-F0kT * std::pow(1.0 -
                                                   std::pow(std::fabs(tau/strength),A_),B_))
                       / std::pow(strength, 3.0) * std::copysign(1.0, tau);
}

}
