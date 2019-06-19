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
  history.add_scalar("strength");
}

void SlipSingleStrengthHardening::init_history(History & history) const
{
  history.get_scalar("strength") = init_strength();
}

History SlipSingleStrengthHardening::hist(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    const Lattice & L, double T, const SlipRule & R) const
{
  double strength = history.get_scalar("strength");

  History res;
  res.add_scalar("strength");
  res.get_scalar("strength") = hist_rate(stress, Q, strength, history, L, T, R);

  return res;
}

History SlipSingleStrengthHardening::d_hist_d_s(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    const Lattice & L, double T,
    const SlipRule & R) const
{
  double strength = history.get_scalar("strength");

  History res;
  res.add_object<Symmetric>("strength");
  res.get_object<Symmetric>("strength") = d_hist_rate_d_stress(stress,
                                                               Q, strength,
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
  double strength = history.get_scalar("strength");

  History res;
  res.add_scalar("strength");
  res.get_scalar("strength") = d_hist_rate_d_strength(stress, Q, strength, 
                                                      history, L, T, R);
  return res;
}

double SlipSingleStrengthHardening::hist_map(const History & history) const
{
  return history.get_scalar("strength");
}

History SlipSingleStrengthHardening::d_hist_map(const History & history) const
{
  History res;
  res.add_scalar("strength");
  res.get_scalar("strength") = 1.0;
  return res;
}

double PlasticSlipHardening::hist_rate(
    const Symmetric & stress, const Orientation & Q, double strength,
    const History & history, const Lattice & L, double T, const SlipRule & R) const
{
  return hist_factor(strength, L, T) * sum_slip_(stress, Q, history, L, T, R); 
}

Symmetric PlasticSlipHardening::d_hist_rate_d_stress(
    const Symmetric & stress, const Orientation & Q, 
    double strength, const History & history, const Lattice & L, double T,
    const SlipRule & R) const
{

}

double PlasticSlipHardening::d_hist_rate_d_strength(
    const Symmetric & stress, const Orientation & Q, 
    double strength, const History & history, const Lattice & L, double T,
    const SlipRule & R) const
{

}

double PlasticSlipHardening::sum_slip_(
    const Symmetric & stress, const Orientation & Q, const History & history,
    const Lattice & L, double T, const SlipRule & R) const
{
  double dg = 0.0;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      dg += R.slip(g, i, stress, Q, history, L, T);
    }
  }

  return dg;
}

Symmetric PlasticSlipHardening::d_sum_slip_d_stress_(
    const Symmetric & stress, const Orientation & Q, const History & history,
    const Lattice & L, double T, const SlipRule & R) const
{
  Symmetric ds;
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      g = R.slip(g, i, stress, Q, history, L, T);
      ds += copysign(1.0, g) * R.d_slip_d_s(g, i, stress, Q, history, L, T);
    }
  }

  return dg;
}

} // namespace neml
