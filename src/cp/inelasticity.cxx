#include "inelasticity.h"

namespace neml {

NoInelasticity::NoInelasticity()
{

}

NoInelasticity::~NoInelasticity()
{

}

std::string NoInelasticity::type()
{
  return "NoInelasticity";
}

ParameterSet NoInelasticity::parameters()
{
  ParameterSet pset(NoInelasticity::type());

  return pset;
}

std::unique_ptr<NEMLObject> NoInelasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<NoInelasticity>(); 
}

void NoInelasticity::populate_history(History & history) const
{
  return;
}

void NoInelasticity::init_history(History & history) const
{
  return;
}

Symmetric NoInelasticity::d_p(const Symmetric & stress, const Orientation & Q,
                              const History & history,
                              const Lattice & lattice, double T) const
{
  return Symmetric();
}

SymSym NoInelasticity::d_d_p_d_stress(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    const Lattice & lattice, double T) const
{
  return SymSym();
}

History NoInelasticity::d_d_p_d_history(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    const Lattice & lattice, double T) const
{
  return History();
}

History NoInelasticity::history_rate(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history,
                                     const Lattice & lattice, double T) const
{
  return History();
}

History NoInelasticity::d_history_rate_d_stress(const Symmetric & stress, 
                                                const Orientation & Q,
                                                const History & history,
                                                const Lattice & lattice, double T) const
{
  return History();
}

History NoInelasticity::d_history_rate_d_history(const Symmetric & stress,
                                               const Orientation & Q,
                                               const History & history,
                                               const Lattice & lattice, double T) const
{
  return History();
}

Skew NoInelasticity::w_p(const Symmetric & stress, const Orientation & Q,
                         const History & history,
                         const Lattice & lattice, double T) const
{
  return Skew();
}


AsaroInelasticity::AsaroInelasticity(std::shared_ptr<SlipRule> rule) :
    rule_(rule)
{

}

AsaroInelasticity::~AsaroInelasticity()
{

}

std::string AsaroInelasticity::type()
{
  return "AsaroInelasticity";
}

ParameterSet AsaroInelasticity::parameters()
{
  ParameterSet pset(AsaroInelasticity::type());

  pset.add_parameter<NEMLObject>("rule");

  return pset;
}

std::unique_ptr<NEMLObject> AsaroInelasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<AsaroInelasticity>(
      params.get_object_parameter<SlipRule>("rule")
      ); 
}

void AsaroInelasticity::populate_history(History & history) const
{
  rule_->populate_history(history);
}

void AsaroInelasticity::init_history(History & history) const
{
  rule_->init_history(history);
}

Symmetric AsaroInelasticity::d_p(const Symmetric & stress, const Orientation & Q,
                              const History & history,
                              const Lattice & lattice, double T) const
{
  Symmetric d;
  for (size_t g = 0; g < lattice.ngroup(); g++) {
    for (size_t i = 0; i < lattice.nslip(g); i++) {
      d += rule_->slip(g, i, stress, Q, history, lattice, T) * lattice.M(
          g, i, Q);
    }
  }

  return d;
}

SymSym AsaroInelasticity::d_d_p_d_stress(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    const Lattice & lattice, double T) const
{
  SymSym ds;

  for (size_t g = 0; g < lattice.ngroup(); g++) {
    for (size_t i = 0; i < lattice.nslip(g); i++) {
      ds += douter(rule_->d_slip_d_s(g, i, stress, Q, history, lattice, T),
                   lattice.M(g, i, Q));
    }
  }

  return ds;
}

History AsaroInelasticity::d_d_p_d_history(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    const Lattice & lattice, double T) const
{
  History h = history.derivative<double>().derivative<Symmetric>();

  for (size_t g = 0; g < lattice.ngroup(); g++) {
    for (size_t i = 0; i < lattice.nslip(g); i++) {
      History h_gi = rule_->d_slip_d_h(g, i, stress, Q, history, lattice, T);
      for (auto item : h_gi.items()) {
        // God is this annoying
        h.get<Symmetric>(item) += h_gi.get<double>(item) * lattice.M(g, i, Q);
      }
    }
  }

  return h;
}

History AsaroInelasticity::history_rate(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history,
                                     const Lattice & lattice, double T) const
{
  return rule_->hist_rate(stress, Q, history, lattice, T);
}

History AsaroInelasticity::d_history_rate_d_stress(const Symmetric & stress, 
                                                const Orientation & Q,
                                                const History & history,
                                                const Lattice & lattice, double T) const
{
  return rule_->d_hist_rate_d_stress(stress, Q, history, lattice, T);
}

History AsaroInelasticity::d_history_rate_d_history(const Symmetric & stress,
                                               const Orientation & Q,
                                               const History & history,
                                               const Lattice & lattice, double T) const
{
  return rule_->d_hist_rate_d_hist(stress, Q, history, lattice, T);
}

Skew AsaroInelasticity::w_p(const Symmetric & stress, const Orientation & Q,
                         const History & history,
                         const Lattice & lattice, double T) const
{
  Skew w;
  for (size_t g = 0; g < lattice.ngroup(); g++) {
    for (size_t i = 0; i < lattice.nslip(g); i++) {
      w += rule_->slip(g, i, stress, Q, history, lattice, T) * lattice.N(
          g, i, Q);
    }
  }

  return w;
}

} // namespace neml
