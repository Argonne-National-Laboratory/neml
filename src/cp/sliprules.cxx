#include "sliprules.h"

namespace neml {

double SlipRule::sum_slip(const Symmetric & stress, const Orientation & Q, 
                          const History & history, const Lattice & L, 
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
                                        const Lattice & L, double T) const
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
                                    const History & history, const Lattice & L,
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

SlipStrengthSlipRule::SlipStrengthSlipRule(
    std::shared_ptr<SlipHardening> strength) :
      strength_(strength)
{
  
}

void SlipStrengthSlipRule::populate_history(History & history) const
{
  strength_->populate_history(history);
}

void SlipStrengthSlipRule::init_history(History & history) const
{
  strength_->init_history(history);
}

double SlipStrengthSlipRule::slip(size_t g, size_t i, const Symmetric & stress, 
                    const Orientation & Q, const History & history,
                    const Lattice & L, double T) const
{
  double tau = L.shear(g, i, Q, stress);
  double tau_bar = strength_->hist_to_tau(g, i, history);

  return sslip(g, i, tau, tau_bar, T);
}

Symmetric SlipStrengthSlipRule::d_slip_d_s(size_t g, size_t i, const Symmetric & stress, 
                             const Orientation & Q, const History & history,
                             const Lattice & L, double T) const
{
  double tau = L.shear(g, i, Q, stress);  
  Symmetric dtau = L.d_shear(g, i, Q, stress);
  double tau_bar = strength_->hist_to_tau(g, i, history);

  return d_sslip_dtau(g, i, tau, tau_bar, T) * dtau;
}

History SlipStrengthSlipRule::d_slip_d_h(size_t g, size_t i, const Symmetric & stress, 
                   const Orientation & Q, const History & history,
                   const Lattice & L, double T) const
{
  double tau = L.shear(g, i, Q, stress);  
  Symmetric dtau = L.d_shear(g, i, Q, stress);
  double tau_bar = strength_->hist_to_tau(g, i, history);

  double dtb = d_sslip_dstrength(g, i, tau, tau_bar, T);

  History deriv = strength_->d_hist_to_tau(g, i, history);
  deriv.scalar_multiply(dtb);

  return deriv;
}

History SlipStrengthSlipRule::hist_rate(const Symmetric & stress, 
                    const Orientation & Q, const History & history,
                    const Lattice & L, double T) const
{
  return strength_->hist(stress, Q, history, L, T, *this);
}

History SlipStrengthSlipRule::d_hist_rate_d_stress(const Symmetric & stress, 
                    const Orientation & Q, const History & history,
                    const Lattice & L, double T) const
{
  return strength_->d_hist_d_s(stress, Q, history, L, T, *this);
}

History SlipStrengthSlipRule::d_hist_rate_d_hist(const Symmetric & stress, 
                    const Orientation & Q, const History & history,
                    const Lattice & L, double T) const
{
  return strength_->d_hist_d_h(stress, Q, history, L, T, *this);
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

double PowerLawSlipRule::sslip(size_t g, size_t i, double tau, double strength, 
                               double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  return g0 * tau / strength * pow(fabs(tau/strength), n-1.0);
}

double PowerLawSlipRule::d_sslip_dtau(size_t g, size_t i, double tau, double strength, 
                                      double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  return g0 * n * pow(fabs(tau/strength), n-1.0) / strength;
}

double PowerLawSlipRule::d_sslip_dstrength(size_t g, size_t i, double tau,
                                           double strength, double T) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);

  return -n * g0 * tau * pow(fabs(tau), n -1.0) / pow(strength, n + 1.0); 
}

} // namespace neml
