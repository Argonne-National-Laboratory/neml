#include "sliprules.h"

namespace neml {

double SlipRule::sum_slip(const Symmetric & stress, const Orientation & Q,
                          const History & history, Lattice & L,
                          double T, const History & fixed) const
{
  double dg = 0.0;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      dg += std::fabs(slip(g, i, stress, Q, history, L, T, fixed));
    }
  }

  return dg;
}

Symmetric SlipRule::d_sum_slip_d_stress(const Symmetric & stress,
                                        const Orientation & Q,
                                        const History & history,
                                        Lattice & L, double T,
                                        const History & fixed) const
{
  Symmetric ds;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      double dg = slip(g, i, stress, Q, history, L, T, fixed);
      ds += std::copysign(1.0, dg) * d_slip_d_s(g, i, stress, Q, history, L, T,
                                           fixed);
    }
  }

  return ds;
}

History SlipRule::d_sum_slip_d_hist(const Symmetric & stress,
                                    const Orientation & Q,
                                    const History & history, Lattice & L,
                                    double T, const History & fixed) const
{
  History res = history.copy_blank();
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      double dg = slip(g, i, stress, Q, history, L, T, fixed);
      double sgn = std::copysign(1.0, dg);
      History temp = d_slip_d_h(g, i, stress, Q,  history, L, T, fixed);
      temp.scalar_multiply(sgn);
      res += temp;
    }
  }

  return res;
}

bool SlipRule::use_nye() const
{
  return false;
}

SlipMultiStrengthSlipRule::SlipMultiStrengthSlipRule(
    std::vector<std::shared_ptr<SlipHardening>> strengths) :
      strengths_(strengths)
{
  // Rename variables to avoid conflicts
  if (strengths.size() > 1) {
    for (size_t i = 0; i < strengths_.size(); i++) {
      auto vars = strengths_[i]->varnames();
      for (size_t j = 0; j < vars.size(); j++) {
        vars[j] += ("_#"+std::to_string(i));
      }
      strengths_[i]->set_varnames(vars);
    }
  }
}

size_t SlipMultiStrengthSlipRule::nstrength() const
{
  return strengths_.size();
}

void SlipMultiStrengthSlipRule::populate_history(History & history) const
{
  for (auto strength : strengths_) {
    strength->populate_history(history);
  }
}

void SlipMultiStrengthSlipRule::init_history(History & history) const
{
  for (auto strength : strengths_) {
    strength->init_history(history);
  }
}

double SlipMultiStrengthSlipRule::strength(const History & history,
                                      Lattice & L, double T,
                                      const History & fixed) const
{
  double val = 0.0;
  double num = 0.0;
  for (auto strength : strengths_) {
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t i = 0; i < L.nslip(g); i++) {
        val += strength->hist_to_tau(g, i, history, L, T, fixed);
        num += 1.0;
      }
    }
  }

  return val / num;
}

double SlipMultiStrengthSlipRule::slip(size_t g, size_t i, const Symmetric & stress,
                    const Orientation & Q, const History & history,
                    Lattice & L, double T, const History & fixed) const
{
  double tau = L.shear(g, i, Q, stress);
  std::vector<double> strengths(nstrength());

  for (size_t i = 0; i < nstrength(); i++) {
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, L, T, fixed);
  }

  return sslip(g, i, tau, strengths, T);
}

Symmetric SlipMultiStrengthSlipRule::d_slip_d_s(size_t g, size_t i, const Symmetric & stress,
                             const Orientation & Q, const History & history,
                             Lattice & L, double T, const History & fixed) const
{
  double tau = L.shear(g, i, Q, stress);
  Symmetric dtau = L.d_shear(g, i, Q, stress);

  std::vector<double> strengths(nstrength());
  for (size_t i = 0; i < nstrength(); i++) {
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, L, T, fixed);
  }

  return d_sslip_dtau(g, i, tau, strengths, T) * dtau;
}

History SlipMultiStrengthSlipRule::d_slip_d_h(size_t g, size_t i, const Symmetric & stress,
                   const Orientation & Q, const History & history,
                   Lattice & L, double T, const History & fixed) const
{
  double tau = L.shear(g, i, Q, stress);
  Symmetric dtau = L.d_shear(g, i, Q, stress);

  std::vector<double> strengths(nstrength());
  for (size_t i = 0; i < nstrength(); i++) {
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, L, T, fixed);
  }

  std::vector<double> dtb = d_sslip_dstrength(g, i, tau, strengths, T);

  History dhist;
  for (size_t i = 0; i < nstrength(); i++) {
    History deriv = strengths_[i]->d_hist_to_tau(g, i, history, L, T, fixed);
    deriv.scalar_multiply(dtb[i]);
    dhist.add_union(deriv);
  }

  return dhist;
}

History SlipMultiStrengthSlipRule::hist_rate(const Symmetric & stress,
                    const Orientation & Q, const History & history,
                    Lattice & L, double T, const History & fixed) const
{
  History res;
  for (auto strength : strengths_) {
    res.add_union(strength->hist(stress, Q, history, L, T, *this, fixed));
  }
  return res;
}

History SlipMultiStrengthSlipRule::d_hist_rate_d_stress(const Symmetric & stress,
                    const Orientation & Q, const History & history,
                    Lattice & L, double T, const History & fixed) const
{
  History dhist;
  for (auto strength : strengths_) {
    dhist.add_union(strength->d_hist_d_s(stress, Q, history, L, T, *this, fixed));
  }
  return dhist;
}

History SlipMultiStrengthSlipRule::d_hist_rate_d_hist(const Symmetric & stress,
                    const Orientation & Q, const History & history,
                    Lattice & L, double T, const History & fixed) const
{
  std::vector<std::string> names;
  History dhist;
  for (size_t i = 0; i < nstrength(); i++) {
    auto vns = strengths_[i]->varnames();
    names.insert(names.end(), vns.begin(), vns.end());
    for (size_t j = 0; j < nstrength(); j++) {
      if (i == j) {
        dhist.add_union(strengths_[i]->d_hist_d_h(stress, Q, history, L, T, *this, fixed));
      }
      else {
        dhist.add_union(strengths_[i]->d_hist_d_h_ext(stress, Q, history, L, T,
                                                      *this, fixed,
                                                      strengths_[j]->varnames()));
      }
    }
  }

  // Most regrettably this is now entirely out of order.  Need to fix up so
  // things are row-major
  std::vector<std::string> order;
  for (auto n1 :  names) {
    for (auto n2 : names) {
      order.push_back(n1 + "_" + n2);
    }
  }
  dhist.reorder(order);

  return dhist;
}

bool SlipMultiStrengthSlipRule::use_nye() const
{
  for (auto strength : strengths_) {
    if (strength->use_nye()) return true;
  }
  return false;
}

KinematicPowerLawSlipRule::KinematicPowerLawSlipRule(
    std::shared_ptr<SlipHardening> backstrength,
    std::shared_ptr<SlipHardening> isostrength,
    std::shared_ptr<SlipHardening> flowresistance,
    std::shared_ptr<Interpolate> gamma0,
    std::shared_ptr<Interpolate> n) :
      SlipMultiStrengthSlipRule({backstrength, isostrength, flowresistance}),
      gamma0_(gamma0), n_(n)
{

}

std::string KinematicPowerLawSlipRule::type()
{
  return "KinematicPowerLawSlipRule";
}

std::unique_ptr<NEMLObject> KinematicPowerLawSlipRule::initialize(
    ParameterSet & params)
{
  return neml::make_unique<KinematicPowerLawSlipRule>(
      params.get_object_parameter<SlipHardening>("backstrength"),
      params.get_object_parameter<SlipHardening>("isostrength"),
      params.get_object_parameter<SlipHardening>("flowresistance"),
      params.get_object_parameter<Interpolate>("gamma0"),
      params.get_object_parameter<Interpolate>("n"));
}

ParameterSet KinematicPowerLawSlipRule::parameters()
{
  ParameterSet pset(KinematicPowerLawSlipRule::type());

  pset.add_parameter<NEMLObject>("backstrength");
  pset.add_parameter<NEMLObject>("isostrength");
  pset.add_parameter<NEMLObject>("flowresistance");
  pset.add_parameter<NEMLObject>("gamma0");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

double KinematicPowerLawSlipRule::sslip(size_t g, size_t i, double tau,
                                        std::vector<double> strengths,
                                        double T) const
{
  double bs = strengths[0];
  double is = strengths[1];
  double fr = strengths[2];

  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  if ((std::fabs(tau - bs) - is) <= 0.0) {
    return 0.0;
  }
  else {
    return std::copysign(g0 * std::pow((std::fabs(tau - bs) - is) / fr, n), tau-bs);
  }
}

double KinematicPowerLawSlipRule::d_sslip_dtau(size_t g, size_t i, double tau,
                                               std::vector<double> strengths,
                                               double T) const
{
  double bs = strengths[0];
  double is = strengths[1];
  double fr = strengths[2];

  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  if ((std::fabs(tau - bs) - is) <= 0.0) {
    return 0.0;
  }
  else {
    return g0 * n * std::pow((std::fabs(tau - bs) - is) / fr, n-1) / fr;
  }
}

std::vector<double> KinematicPowerLawSlipRule::d_sslip_dstrength(
    size_t g, size_t i, double tau, std::vector<double> strengths,
    double T) const
{
  double bs = strengths[0];
  double is = strengths[1];
  double fr = strengths[2];

  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  if ((std::fabs(tau - bs) - is) <= 0.0) {
    return {0.0,0.0,0.0};
  }
  else {
    double dbs = g0 * n * std::pow((std::fabs(tau-bs) - is) / fr, n) / (is - std::fabs(tau - bs));
    double dis = -std::copysign(g0 * n * std::pow((std::fabs(tau-bs) - is) / fr, n-1) / fr, tau - bs);
    double dfr = -std::copysign(g0 * n * std::pow((std::fabs(tau-bs) - is) / fr, n) / fr, tau - bs);
    return {dbs, dis, dfr};
  }
}

SlipStrengthSlipRule::SlipStrengthSlipRule(
    std::shared_ptr<SlipHardening> strength) :
      SlipMultiStrengthSlipRule({strength})
{

}

double SlipStrengthSlipRule::sslip(size_t g, size_t i, double tau,
                                   std::vector<double> strengths, double T) const
{
  return scalar_sslip(g, i, tau, strengths[0], T);
}

double SlipStrengthSlipRule::d_sslip_dtau(size_t g, size_t i, double tau,
                                          std::vector<double> strengths,
                                          double T) const
{
  return scalar_d_sslip_dtau(g, i, tau, strengths[0], T);
}

std::vector<double> SlipStrengthSlipRule::d_sslip_dstrength(size_t g,
                                                            size_t i,
                                                            double tau,
                                                            std::vector<double> strengths,
                                                            double T) const
{
  return {scalar_d_sslip_dstrength(g, i, tau, strengths[0], T)};
}


HuCocksSlipRule::HuCocksSlipRule(std::shared_ptr<SlipHardening> strength,
                                   std::shared_ptr<Interpolate> gamma0,
                                   std::shared_ptr<Interpolate> n) :
    SlipStrengthSlipRule(strength), gamma0_(gamma0), n_(n)
{

}

std::string HuCocksSlipRule::type()
{
  return "HuCocksSlipRule";
}

std::unique_ptr<NEMLObject> HuCocksSlipRule::initialize(
    ParameterSet & params)
{
  return neml::make_unique<HuCocksSlipRule>(
      params.get_object_parameter<SlipHardening>("resistance"),
      params.get_object_parameter<Interpolate>("gamma0"),
      params.get_object_parameter<Interpolate>("n"));
}

ParameterSet HuCocksSlipRule::parameters()
{
  ParameterSet pset(HuCocksSlipRule::type());

  pset.add_parameter<NEMLObject>("resistance");
  pset.add_parameter<NEMLObject>("gamma0");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

double HuCocksSlipRule::scalar_sslip(size_t g, size_t i, double tau,
                                      double strength, double T) const
{

  double g0 = gamma0_->value(T);
  double n = n_->value(T);
  double Fo = 0.5 * 15.625 * std::pow(10,-30) * 79.3 * std::pow(10,9);
  double k = 1.38 * std::pow(10,-23);
  double Temp = 1200;
  double C = Fo/k/Temp;
  //tau = 1;
  double tau_in = 1;
 // return g0 * tau / strength * std::pow(std::fabs(tau/strength), n-1.0);
  return g0 * std::exp(-C * std::pow(1 - std::pow((tau + tau_in) / strength,1),1)) * (tau+tau_in)/std::abs(tau+tau_in);
//  return g0 * (tau+tau_in) / strength * std::pow(std::fabs(tau/strength), n-1.0);
//  return 0.01*std::exp(-C * std::pow(1 - std::pow(std::abs((tau + tau_in)/strength),1),1))
//  * (tau/std::abs(tau));

}

double HuCocksSlipRule::scalar_d_sslip_dtau(size_t g, size_t i, double tau,
                                             double strength, double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);
  double Fo = 0.5 * 15.625 * std::pow(10,-30) * 79.3 * std::pow(10,9);
  double k = 1.38 * std::pow(10,-23);
  double Temp = 1200;
  double C = Fo/k/Temp;
  //tau = 1;
  double tau_in = 1;
//  return C * g0 * std::pow(tau/strength,0.75) * std::pow(1 - std::pow(tau/strength,0.75),0.33)
//  * std::exp(-C * std::pow(1 - std::pow(tau/strength,0.75),1.33))/tau;
  return C * g0 * std::exp(-C * (1 - (tau + tau_in)/strength))/strength  * (tau+tau_in)/std::abs(tau+tau_in);
  // return g0 * n * std::pow(std::fabs(tau/strength), n-1.0) / strength;
//  return 0;
//  return g0*C*std::pow((tau + tau_in)/strength,0.75)*std::pow(1 - std::pow((tau + tau_in)/strength,0.75),
//  0.33)*std::exp(-C*std::pow(1 - std::pow((tau + tau_in)/strength,0.75),1.33))/(tau + tau_in);

}

double HuCocksSlipRule::scalar_d_sslip_dstrength(size_t g, size_t i,
                                                  double tau,
                                                  double strength,
                                                  double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);
  double Fo = 0.65 * 15.625 * std::pow(10,-30) * 79.3 * std::pow(10,9);
  double k = 1.38 * std::pow(10,-23);
  double Temp = 1200;
  double C = Fo/k/Temp;
  //tau = 1;
  double tau_in = 1;

//  return 0;
//  return -C * g0 * std::pow(tau/strength,1) * std::pow(1 - std::pow(tau/strength,1),1) *
//  std::exp(-C * std::pow(1 - std::pow(tau/strength,1), 1))/strength;
//  return -C * g0 * std::pow(tau/strength,0.75) * std::pow(1 - std::pow(tau/strength,0.75),0.33)
//  * std::exp(-C * std::pow(1 - std::pow(tau/strength,0.75),1.33))/strength;
  return -C * g0 * (tau + tau_in) * std::exp(-C * (1 - (tau + tau_in)/strength))/ std::pow(strength,2) * (tau+tau_in)/std::abs(tau+tau_in);
  // return -n * g0 * tau * std::pow(std::fabs(tau), n -1.0) / std::pow(strength, n + 1.0);
//  return -g0*C*std::pow((tau + tau_in)/strength,0.75)*std::pow(1 - std::pow((tau + tau_in)/strength,0.75),
//  0.33)*std::exp(-C*std::pow(1 - std::pow((tau + tau_in)/strength,0.75),1.33))/strength;

}


PowerLawSlipRule::PowerLawSlipRule(std::shared_ptr<SlipHardening> strength,
                                   std::shared_ptr<Interpolate> gamma0,
                                   std::shared_ptr<Interpolate> n) :
    SlipStrengthSlipRule(strength), gamma0_(gamma0), n_(n)
{

}

std::string PowerLawSlipRule::type()
{
  return "PowerLawSlipRule";
}

std::unique_ptr<NEMLObject> PowerLawSlipRule::initialize(
    ParameterSet & params)
{
  return neml::make_unique<PowerLawSlipRule>(
      params.get_object_parameter<SlipHardening>("resistance"),
      params.get_object_parameter<Interpolate>("gamma0"),
      params.get_object_parameter<Interpolate>("n"));
}

ParameterSet PowerLawSlipRule::parameters()
{
  ParameterSet pset(PowerLawSlipRule::type());

  pset.add_parameter<NEMLObject>("resistance");
  pset.add_parameter<NEMLObject>("gamma0");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

double PowerLawSlipRule::scalar_sslip(size_t g, size_t i, double tau,
                                      double strength, double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  return g0 * tau / strength * std::pow(std::fabs(tau/strength), n-1.0);
}

double PowerLawSlipRule::scalar_d_sslip_dtau(size_t g, size_t i, double tau,
                                             double strength, double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  return g0 * n * std::pow(std::fabs(tau/strength), n-1.0) / strength;
}

double PowerLawSlipRule::scalar_d_sslip_dstrength(size_t g, size_t i,
                                                  double tau,
                                                  double strength,
                                                  double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  return -n * g0 * tau * std::pow(std::fabs(tau), n -1.0) / std::pow(strength, n + 1.0);
}

} // namespace neml
