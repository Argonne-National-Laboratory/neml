#include "slipharden.h"

namespace neml {

void SingleSlipHardening::populate_history(History & history) const
{
  history.add_scalar("strength");
}

void SingleSlipHardening::init_history(History & history) const
{
  history.get_scalar("strength") = 0.0;
}

double SingleSlipHardening::hist_to_tau(size_t g, size_t i,
                                        const History & history) const
{
  return history.get_scalar("strength");
}

std::map<std::string,std::vector<double>> 
SingleSlipHardening::d_hist_to_tau(size_t g, size_t i,
                                   const History & history) const
{
  std::map<std::string,std::vector<double>> vals;
  vals["strength"] = std::vector<double>({1.0});
  return vals;
}

History SingleSSSlipHardening::hist(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    const Lattice & L, double T, const SlipRule & R) const
{

}

Symmetric d_hist_d_s(const Symmetric & stress, 
                             const Orientation & Q, const History & history,
                             const Lattice & L, double T,
                             const SlipRule & R) const;
std::map<std::string,std::vector<double>>
    d_hist_d_h(const Symmetric & stress, 
               const Orientation & Q,
               const History & history,
               const Lattice & L, 
               double T, const SlipRule & R) const;

} // namespace neml
