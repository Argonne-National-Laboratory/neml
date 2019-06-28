#include "slipharden.h"

namespace neml {

double SlipSingleHardening::hist_to_tau(
    size_t g, size_t i, const History & history) const
{
  return hist_map(history);
}

History SlipSingleHardening::d_hist_to_tau(
    size_t g, size_t i, const History & history) const
{
  return d_hist_map(history);
}

void SlipSingleStrengthHardening::populate_history(History & history) const
{
  history.add_object<double>("strength");
}

void SlipSingleStrengthHardening::init_history(History & history) const
{
  history.get_object<double>("strength") = init_strength();
}

History SlipSingleStrengthHardening::hist(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    const Lattice & L, double T, const SlipRule & R) const
{
  History res;
  res.add_object<double>("strength");
  res.get_object<double>("strength") = hist_rate(stress, Q, history, L, T, R);

  return res;
}

History SlipSingleStrengthHardening::d_hist_d_s(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    const Lattice & L, double T,
    const SlipRule & R) const
{
  History res;
  res.add_object<Symmetric>("strength");
  res.get_object<Symmetric>("strength") = d_hist_rate_d_stress(stress,
                                                               Q, 
                                                               history,L, T, R);
  return res;
}

History SlipSingleStrengthHardening::d_hist_d_h(
    const Symmetric & stress, 
    const Orientation & Q,
    const History & history,
    const Lattice & L, 
    double T, const SlipRule & R) const
{
  History res;
  res.add_object<double>("strength_strength");
  res.get_object<double>("strength_strength") = d_hist_rate_d_strength(
      stress, Q, history, L, T, R);
  return res;
}

double SlipSingleStrengthHardening::hist_map(const History & history) const
{
  return history.get_object<double>("strength");
}

History SlipSingleStrengthHardening::d_hist_map(const History & history) const
{
  History res;
  res.add_object<double>("strength");
  res.get_object<double>("strength") = 1.0;
  return res;
}

double PlasticSlipHardening::hist_rate(
    const Symmetric & stress, const Orientation & Q,
    const History & history, const Lattice & L, double T, const SlipRule & R) const
{
  double strength = history.get_object<double>("strength");
  return hist_factor(strength, L, T) * R.sum_slip(stress, Q, history, L, T);
}

Symmetric PlasticSlipHardening::d_hist_rate_d_stress(
    const Symmetric & stress, const Orientation & Q, 
    const History & history, const Lattice & L, double T,
    const SlipRule & R) const
{
  double strength = history.get_object<double>("strength");
  return hist_factor(strength, L, T) * R.d_sum_slip_d_stress(stress, Q,
                                                            history, L, T);
}

double PlasticSlipHardening::d_hist_rate_d_strength(
    const Symmetric & stress, const Orientation & Q, 
    const History & history, const Lattice & L, double T,
    const SlipRule & R) const
{
  double strength = history.get_object<double>("strength");
  History dhist = R.d_sum_slip_d_hist(stress, Q, history, L, T);
  double dstrength = dhist.get_object<double>("strength");
  return d_hist_factor(strength, L, T) * R.sum_slip(stress, Q, history, L, T) 
      + hist_factor(strength, L, T) * dstrength;
}

VoceSlipHardening::VoceSlipHardening(std::shared_ptr<Interpolate> tau_sat,
                                     std::shared_ptr<Interpolate> b,
                                     std::shared_ptr<Interpolate> tau_0) :
    tau_sat_(tau_sat), b_(b), tau_0_(tau_0)
{
  
}

std::string VoceSlipHardening::type()
{
  return "VoceSlipHardening";
}

std::unique_ptr<NEMLObject> VoceSlipHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<VoceSlipHardening>(
      params.get_object_parameter<Interpolate>("tau_sat"),
      params.get_object_parameter<Interpolate>("b"),
      params.get_object_parameter<Interpolate>("tau_0"));
}

ParameterSet VoceSlipHardening::parameters()
{
  ParameterSet pset(VoceSlipHardening::type());
  
  pset.add_parameter<NEMLObject>("tau_sat");
  pset.add_parameter<NEMLObject>("b");
  pset.add_parameter<NEMLObject>("tau_0");

  return pset;
}

double VoceSlipHardening::init_strength() const
{
  return 0.0;
}

double VoceSlipHardening::hist_factor(double strength, const Lattice & L, double T) const
{
  double tau_sat = tau_sat_->value(T);
  double b = b_->value(T);
  double tau_0 = tau_0_->value(T);

  return b * (tau_sat - strength) + tau_0;
}

double VoceSlipHardening::d_hist_factor(double strength, const Lattice & L, double T) const
{
  double b = b_->value(T);

  return -b;
}

} // namespace neml
