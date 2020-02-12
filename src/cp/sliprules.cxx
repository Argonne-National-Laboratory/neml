#include "sliprules.h"

namespace neml {

double SlipRule::sum_slip(const Symmetric & stress, const Orientation & Q, 
                          const History & history, Lattice & L, 
                          double T) const
{
  double dg = 0.0;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      dg += fabs(slip(g, i, stress, Q, history, L, T));
    }
  }

  return dg;
}

Symmetric SlipRule::d_sum_slip_d_stress(const Symmetric & stress, 
                                        const Orientation & Q, 
                                        const History & history,
                                        Lattice & L, double T) const
{
  Symmetric ds;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      double dg = slip(g, i, stress, Q, history, L, T);
      ds += copysign(1.0, dg) * d_slip_d_s(g, i, stress, Q, history, L, T);
    }
  }

  return ds;
}

History SlipRule::d_sum_slip_d_hist(const Symmetric & stress,
                                    const Orientation & Q, 
                                    const History & history, Lattice & L,
                                    double T) const
{
  History res = history.copy_blank();
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      double dg = slip(g, i, stress, Q, history, L, T);
      double sgn = copysign(1.0, dg);
      History temp = d_slip_d_h(g, i, stress, Q,  history, L, T);
      temp.scalar_multiply(sgn);
      res += temp;
    }
  }

  return res;
}

SlipMultiStrengthSlipRule::SlipMultiStrengthSlipRule(
    std::vector<std::shared_ptr<SlipHardening>> strengths) :
      strengths_(strengths)
{
  
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
                                      Lattice & L, double T) const
{
  double val = 0.0;
  double num = 0.0;
  for (auto strength : strengths_) {
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t i = 0; i < L.nslip(g); i++) {
        val += strength->hist_to_tau(g, i, history, T);
        num += 1.0;
      }
    }
  }

  return val / num;
}

double SlipMultiStrengthSlipRule::slip(size_t g, size_t i, const Symmetric & stress, 
                    const Orientation & Q, const History & history,
                    Lattice & L, double T) const
{
  double tau = L.shear(g, i, Q, stress);
  std::vector<double> strengths(nstrength());
  for (size_t i = 0; i < nstrength(); i++) {
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, T);
  }

  return sslip(g, i, tau, strengths, T);
}

Symmetric SlipMultiStrengthSlipRule::d_slip_d_s(size_t g, size_t i, const Symmetric & stress, 
                             const Orientation & Q, const History & history,
                             Lattice & L, double T) const
{
  double tau = L.shear(g, i, Q, stress);  
  Symmetric dtau = L.d_shear(g, i, Q, stress);

  std::vector<double> strengths(nstrength());
  for (size_t i = 0; i < nstrength(); i++) {
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, T);
  }

  return d_sslip_dtau(g, i, tau, strengths, T) * dtau;
}

// This is the hard one
History SlipMultiStrengthSlipRule::d_slip_d_h(size_t g, size_t i, const Symmetric & stress, 
                   const Orientation & Q, const History & history,
                   Lattice & L, double T) const
{
  double tau = L.shear(g, i, Q, stress);  
  Symmetric dtau = L.d_shear(g, i, Q, stress);
  std::vector<double> strengths(nstrength());
  for (size_t i = 0; i < nstrength(); i++) {
    strengths[i] = strengths_[i]->hist_to_tau(g, i, history, T);
  }

  std::vector<double> dtb = d_sslip_dstrength(g, i, tau, strengths, T);

  History dhist;

  for (size_t i = 0; i < nstrength(); i++) {
    History deriv = strengths_[i]->d_hist_to_tau(g, i, history, T);
    deriv.scalar_multiply(dtb[i]);
    dhist.add_union(deriv); 
  }

  return dhist;
}

History SlipMultiStrengthSlipRule::hist_rate(const Symmetric & stress, 
                    const Orientation & Q, const History & history,
                    Lattice & L, double T) const
{
  History hist;
  for (auto strength : strengths_) {
    hist.add_union(strength->hist(stress, Q, history, L, T, *this));
  }
  return hist;
}

History SlipMultiStrengthSlipRule::d_hist_rate_d_stress(const Symmetric & stress, 
                    const Orientation & Q, const History & history,
                    Lattice & L, double T) const
{
  History dhist;
  for (auto strength : strengths_) {
    dhist.add_union(strength->d_hist_d_s(stress, Q, history, L, T, *this));
  }
  return dhist;
}

History SlipMultiStrengthSlipRule::d_hist_rate_d_hist(const Symmetric & stress, 
                    const Orientation & Q, const History & history,
                    Lattice & L, double T) const
{
  History dhist;
  for (auto strength : strengths_) {
    dhist.add_union(strength->d_hist_d_h(stress, Q, history, L, T, *this));
  }
  return dhist;
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
