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
                              Lattice & lattice, double T) const
{
  return Symmetric();
}

SymSym NoInelasticity::d_d_p_d_stress(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  return SymSym();
}

History NoInelasticity::d_d_p_d_history(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  return History();
}

History NoInelasticity::history_rate(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history,
                                     Lattice & lattice, double T) const
{
  return History();
}

History NoInelasticity::d_history_rate_d_stress(const Symmetric & stress, 
                                                const Orientation & Q,
                                                const History & history,
                                                Lattice & lattice, double T) const
{
  return History();
}

History NoInelasticity::d_history_rate_d_history(const Symmetric & stress,
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & lattice, double T) const
{
  return History();
}

Skew NoInelasticity::w_p(const Symmetric & stress, const Orientation & Q,
                         const History & history,
                         Lattice & lattice, double T) const
{
  return Skew();
}

SkewSym NoInelasticity::d_w_p_d_stress(const Symmetric & stress, 
                                       const Orientation & Q,
                                       const History & history,
                                       Lattice & lattice, double T) const
{
  return SkewSym();
}

History NoInelasticity::d_w_p_d_history(const Symmetric & stress,
                                        const Orientation & Q,
                                        const History & history,
                                        Lattice & lattice, double T) const
{
  return History();
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
                              Lattice & lattice, double T) const
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
    Lattice & lattice, double T) const
{
  SymSym ds;

  for (size_t g = 0; g < lattice.ngroup(); g++) {
    for (size_t i = 0; i < lattice.nslip(g); i++) {
      ds += douter(lattice.M(g, i, Q),
                   rule_->d_slip_d_s(g, i, stress, Q, history, lattice, T));
    }
  }

  return ds;
}

History AsaroInelasticity::d_d_p_d_history(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  History h = history.derivative<Symmetric>();

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
                                     Lattice & lattice, double T) const
{
  return rule_->hist_rate(stress, Q, history, lattice, T);
}

History AsaroInelasticity::d_history_rate_d_stress(const Symmetric & stress, 
                                                const Orientation & Q,
                                                const History & history,
                                                Lattice & lattice, double T) const
{
  return rule_->d_hist_rate_d_stress(stress, Q, history, lattice, T);
}

History AsaroInelasticity::d_history_rate_d_history(const Symmetric & stress,
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & lattice, double T) const
{
  return rule_->d_hist_rate_d_hist(stress, Q, history, lattice, T);
}

Skew AsaroInelasticity::w_p(const Symmetric & stress, const Orientation & Q,
                         const History & history,
                         Lattice & lattice, double T) const
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

SkewSym AsaroInelasticity::d_w_p_d_stress(const Symmetric & stress, 
                                       const Orientation & Q,
                                       const History & history,
                                       Lattice & lattice, double T) const
{
  SkewSym ds;

  for (size_t g = 0; g < lattice.ngroup(); g++) {
    for (size_t i = 0; i < lattice.nslip(g); i++) {
      ds += douter(lattice.N(g, i, Q), 
                   rule_->d_slip_d_s(g, i, stress, Q, history, lattice, T));
    }
  }

  return ds;
}

History AsaroInelasticity::d_w_p_d_history(const Symmetric & stress,
                                        const Orientation & Q,
                                        const History & history,
                                        Lattice & lattice, double T) const
{
  History h = history.derivative<Skew>();

  for (size_t g = 0; g < lattice.ngroup(); g++) {
    for (size_t i = 0; i < lattice.nslip(g); i++) {
      History h_gi = rule_->d_slip_d_h(g, i, stress, Q, history, lattice, T);
      for (auto item : h_gi.items()) {
        // God is this annoying
        h.get<Skew>(item) += h_gi.get<double>(item) * lattice.N(g, i, Q);
      }
    }
  }

  return h;
}

PowerLaw::PowerLaw(std::shared_ptr<Interpolate> A, 
                   std::shared_ptr<Interpolate> n) :
    A_(A), n_(n)
{

}

PowerLaw::~PowerLaw()
{

}

std::string PowerLaw::type()
{
  return "PowerLaw";
}

ParameterSet PowerLaw::parameters()
{
  ParameterSet pset(PowerLaw::type());
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

std::unique_ptr<NEMLObject> PowerLaw::initialize(ParameterSet & params)
{
  return neml::make_unique<PowerLaw>(
      params.get_object_parameter<Interpolate>("A"),
      params.get_object_parameter<Interpolate>("n")); 
}

void PowerLaw::populate_history(History & history) const
{
  return;
}

void PowerLaw::init_history(History & history) const
{
  return;
}

Symmetric PowerLaw::d_p(const Symmetric & stress, const Orientation & Q,
                              const History & history,
                              Lattice & lattice, double T) const
{
  double seq = seq_(stress);
  double A = A_->value(T);
  double n = n_->value(T);

  return A * pow(seq, n-1.0) * stress.dev();
}

SymSym PowerLaw::d_d_p_d_stress(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  double seq = seq_(stress);
  double A = A_->value(T);
  double n = n_->value(T);

  Symmetric dir = stress.dev() / seq;
  
  double rate = A * pow(seq, n);
  double drate = A * n * pow(seq, n-1.0);

  return SymSym::id_dev().dot(3.0/2.0 * (drate - rate / seq) * douter(dir,dir) + 
                             SymSym::id() * rate / seq);
}

History PowerLaw::d_d_p_d_history(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  return History();
}

History PowerLaw::history_rate(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history,
                                     Lattice & lattice, double T) const
{
  return History();
}

History PowerLaw::d_history_rate_d_stress(const Symmetric & stress, 
                                                const Orientation & Q,
                                                const History & history,
                                                Lattice & lattice, double T) const
{
  return History();
}

History PowerLaw::d_history_rate_d_history(const Symmetric & stress,
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & lattice, double T) const
{
  return History();
}

Skew PowerLaw::w_p(const Symmetric & stress, const Orientation & Q,
                         const History & history,
                         Lattice & lattice, double T) const
{
  return Skew();
}

SkewSym PowerLaw::d_w_p_d_stress(const Symmetric & stress, 
                                       const Orientation & Q,
                                       const History & history,
                                       Lattice & lattice, double T) const
{
  return SkewSym();
}

History PowerLaw::d_w_p_d_history(const Symmetric & stress,
                                        const Orientation & Q,
                                        const History & history,
                                        Lattice & lattice, double T) const
{
  return History();
}

double PowerLaw::seq_(const Symmetric & stress) const
{
  return sqrt(3.0/2.0) * stress.dev().norm();
}

} // namespace neml
