#include "sliprules.h"

namespace neml {

double SlipRule::sum_slip(const Symmetric & stress, const Orientation & Q, 
                          const History & history, Lattice & L, 
                          double T, const History & fixed) const
{
  double dg = 0.0;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      dg += fabs(slip(g, i, stress, Q, history, L, T, fixed));
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
      ds += copysign(1.0, dg) * d_slip_d_s(g, i, stress, Q, history, L, T,
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
      double sgn = copysign(1.0, dg);
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
        vars[j] += std::to_string(i);
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
        val += strength->hist_to_tau(g, i, history, T, fixed);
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
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, T, fixed);
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
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, T, fixed);
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
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, T, fixed);
  }

  std::vector<double> dtb = d_sslip_dstrength(g, i, tau, strengths, T);

  History dhist;
  for (size_t i = 0; i < nstrength(); i++) {
    History deriv = strengths_[i]->d_hist_to_tau(g, i, history, T, fixed);
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
  History dhist;
  for (auto strength : strengths_) {
    dhist.add_union(strength->d_hist_d_h(stress, Q, history, L, T, *this, fixed));
  }
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
    std::shared_ptr<SlipHardening> understrength,
    std::shared_ptr<Interpolate> gamma0,
    std::shared_ptr<Interpolate> n) :
      SlipMultiStrengthSlipRule({backstrength, understrength}), gamma0_(gamma0),
      n_(n)
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
      params.get_object_parameter<SlipHardening>("understrength"),
      params.get_object_parameter<Interpolate>("gamma0"),
      params.get_object_parameter<Interpolate>("n"));
}

ParameterSet KinematicPowerLawSlipRule::parameters()
{
  ParameterSet pset(KinematicPowerLawSlipRule::type());
  
  pset.add_parameter<NEMLObject>("backstrength");
  pset.add_parameter<NEMLObject>("understrength");
  pset.add_parameter<NEMLObject>("gamma0");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

double KinematicPowerLawSlipRule::sslip(size_t g, size_t i, double tau, 
                                        std::vector<double> strengths, 
                                        double T) const
{
  double bs = strengths[0];
  double us = strengths[1];

  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  if ((tau - bs) < 0.0) {
    return 0.0;
  }
  else {
    return g0 * (tau - bs) / us * pow(fabs((tau-bs)/us), n-1.0);
  }
}

double KinematicPowerLawSlipRule::d_sslip_dtau(size_t g, size_t i, double tau, 
                                               std::vector<double> strengths,
                                               double T) const
{
  double bs = strengths[0];
  double us = strengths[1];

  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  if ((tau - bs) < 0.0) {
    return 0.0;
  }
  else {
    return g0 * n * pow(us,-n) * pow(fabs(bs-tau), n-1.0);
  } 
}

std::vector<double> KinematicPowerLawSlipRule::d_sslip_dstrength(
    size_t g, size_t i, double tau, std::vector<double> strengths, 
    double T) const
{
  double bs = strengths[0];
  double us = strengths[1];

  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  if ((tau - bs) < 0.0) {
    return {0.0,0.0};
  }
  else {
    double d1 = -g0 * n * pow(us,-n) * pow(fabs(bs-tau),n-1.0);
    double d2 = g0 * n * (bs - tau) * pow(us,-1.0-n) * pow(fabs(bs-tau), n - 1.0);
    return {d1, d2};
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

  return g0 * tau / strength * pow(fabs(tau/strength), n-1.0);
}

double PowerLawSlipRule::scalar_d_sslip_dtau(size_t g, size_t i, double tau, 
                                             double strength, double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  return g0 * n * pow(fabs(tau/strength), n-1.0) / strength;
}

double PowerLawSlipRule::scalar_d_sslip_dstrength(size_t g, size_t i, 
                                                  double tau, 
                                                  double strength,
                                                  double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  return -n * g0 * tau * pow(fabs(tau), n -1.0) / pow(strength, n + 1.0); 
}

} // namespace neml
