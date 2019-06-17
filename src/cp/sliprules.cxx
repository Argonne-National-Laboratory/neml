#include "sliprules.h"

namespace neml {

PowerLawSlipRule::PowerLawSlipRule(std::shared_ptr<Interpolate> gamma0, 
                                   std::shared_ptr<Interpolate> n) :
    gamma0_(gamma0), n_(n)
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
      params.get_object_parameter<Interpolate>("gamma0"),
      params.get_object_parameter<Interpolate>("n"));
}

ParameterSet PowerLawSlipRule::parameters()
{
  ParameterSet pset;

  pset.add_parameter<NEMLObject>("gamma0");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

void PowerLawSlipRule::populate_history(History & history) const
{
  return;
}

void PowerLawSlipRule::init_history(History & history) const
{
  return;
}

double PowerLawSlipRule::slip(size_t g, size_t i, const Symmetric & stress, 
                              const Orientation & Q, const History & history,
                              const Lattice & L, double T, 
                              const SlipHardening & H) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);
  double tau = L.shear(g, i, Q, stress);
  double tau_bar = H.hist_to_tau(g, i, history);

  return g0 * tau / tau_bar * pow(fabs(tau/tau_bar), n-1.0);
}

Symmetric PowerLawSlipRule::d_slip_d_s(size_t g, size_t i, 
                                       const Symmetric & stress, 
                                       const Orientation & Q, 
                                       const History & history,
                                       const Lattice & L, 
                                       double T, const SlipHardening & H) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);
  double tau = L.shear(g, i, Q, stress);
  Symmetric dtau = L.d_shear(g, i, Q, stress);
  double tau_bar = H.hist_to_tau(g, i, history);
  
  return g0 * n * pow(fabs(tau/tau_bar), n-1.0) / tau_bar * dtau;

}

std::map<std::string,std::vector<double>> 
PowerLawSlipRule::d_slip_d_h(size_t g, size_t i, 
                             const Symmetric & stress, 
                             const Orientation & Q, 
                             const History & history,
                             const Lattice & L,
                             double T, const SlipHardening & H) const
{
  double g0 = gamma0_->value(T);
  double n = n_->value(T);
  double tau = L.shear(g, i, Q, stress);
  double tau_bar = H.hist_to_tau(g, i, history);

  double sf = -g0 * n * tau * pow(fabs(tau / tau_bar), n-1.0);

  auto lower = H.d_hist_to_tau(g, i, history);
  for (auto it = lower.begin(); it != lower.end(); it++) {
    for (auto sit = it->second.begin(); sit != it->second.end(); sit++) {
      *sit *= sf;
    }
  }

  return lower;
}

} // namespace neml
