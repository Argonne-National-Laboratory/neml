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

double NoInelasticity::strength(const History & history, Lattice & L,
                                double T) const
{
  return 0; // I guess
}

Symmetric NoInelasticity::d_p(const Symmetric & stress, const Orientation & Q,
                              const History & history,
                              Lattice & lattice, double T) const
{
  return Symmetric();
}

SymSymR4 NoInelasticity::d_d_p_d_stress(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  return SymSymR4();
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

SkewSymR4 NoInelasticity::d_w_p_d_stress(const Symmetric & stress, 
                                       const Orientation & Q,
                                       const History & history,
                                       Lattice & lattice, double T) const
{
  return SkewSymR4();
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

double AsaroInelasticity::strength(const History & history, Lattice & L,
                                   double T) const
{
  return rule_->strength(history, L, T);
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

SymSymR4 AsaroInelasticity::d_d_p_d_stress(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  SymSymR4 ds;

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

SkewSymR4 AsaroInelasticity::d_w_p_d_stress(const Symmetric & stress, 
                                       const Orientation & Q,
                                       const History & history,
                                       Lattice & lattice, double T) const
{
  SkewSymR4 ds;

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

PowerLawInelasticity::PowerLawInelasticity(std::shared_ptr<Interpolate> A, 
                   std::shared_ptr<Interpolate> n) :
    A_(A), n_(n)
{

}

PowerLawInelasticity::~PowerLawInelasticity()
{

}

std::string PowerLawInelasticity::type()
{
  return "PowerLawInelasticity";
}

ParameterSet PowerLawInelasticity::parameters()
{
  ParameterSet pset(PowerLawInelasticity::type());
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

std::unique_ptr<NEMLObject> PowerLawInelasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<PowerLawInelasticity>(
      params.get_object_parameter<Interpolate>("A"),
      params.get_object_parameter<Interpolate>("n")); 
}

void PowerLawInelasticity::populate_history(History & history) const
{
  return;
}

void PowerLawInelasticity::init_history(History & history) const
{
  return;
}

double PowerLawInelasticity::strength(const History & history,
                                      Lattice & L, double T) const
{
  return pow(A_->value(T), -1.0/(n_->value(T)));
}

Symmetric PowerLawInelasticity::d_p(const Symmetric & stress, const Orientation & Q,
                              const History & history,
                              Lattice & lattice, double T) const
{
  double seq = seq_(stress);
  double A = A_->value(T);
  double n = n_->value(T);

  if (seq < std::numeric_limits<double>::epsilon()) {
    return Symmetric();
  }

  return A * pow(seq, n-1.0) * stress.dev();
}

SymSymR4 PowerLawInelasticity::d_d_p_d_stress(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  double seq = seq_(stress);
  double A = A_->value(T);
  double n = n_->value(T);

  if (seq < std::numeric_limits<double>::epsilon()) {
    seq = std::numeric_limits<double>::epsilon();
  }

  Symmetric dir = stress.dev() / seq;
  
  double rate = A * pow(seq, n);
  double drate = A * n * pow(seq, n-1.0);

  return SymSymR4::id_dev().dot(3.0/2.0 * (drate - rate / seq) * douter(dir,dir) + 
                             SymSymR4::id() * rate / seq);
}

History PowerLawInelasticity::d_d_p_d_history(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  return History();
}

History PowerLawInelasticity::history_rate(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history,
                                     Lattice & lattice, double T) const
{
  return History();
}

History PowerLawInelasticity::d_history_rate_d_stress(const Symmetric & stress, 
                                                const Orientation & Q,
                                                const History & history,
                                                Lattice & lattice, double T) const
{
  return History();
}

History PowerLawInelasticity::d_history_rate_d_history(const Symmetric & stress,
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & lattice, double T) const
{
  return History();
}

Skew PowerLawInelasticity::w_p(const Symmetric & stress, const Orientation & Q,
                         const History & history,
                         Lattice & lattice, double T) const
{
  return Skew();
}

SkewSymR4 PowerLawInelasticity::d_w_p_d_stress(const Symmetric & stress, 
                                       const Orientation & Q,
                                       const History & history,
                                       Lattice & lattice, double T) const
{
  return SkewSymR4();
}

History PowerLawInelasticity::d_w_p_d_history(const Symmetric & stress,
                                        const Orientation & Q,
                                        const History & history,
                                        Lattice & lattice, double T) const
{
  return History();
}

double PowerLawInelasticity::seq_(const Symmetric & stress) const
{
  return sqrt(3.0/2.0) * stress.dev().norm();
}

CombinedInelasticity::CombinedInelasticity(
    std::vector<std::shared_ptr<InelasticModel>> models) :
      models_(models)
{

}

CombinedInelasticity::~CombinedInelasticity()
{

}

std::string CombinedInelasticity::type()
{
  return "CombinedInelasticity";
}

ParameterSet CombinedInelasticity::parameters()
{
  ParameterSet pset(CombinedInelasticity::type());

  pset.add_parameter<std::vector<NEMLObject>>("models");

  return pset;
}

std::unique_ptr<NEMLObject> CombinedInelasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<CombinedInelasticity>(
      params.get_object_parameter_vector<InelasticModel>("models")
      ); 
}

void CombinedInelasticity::populate_history(History & history) const
{
  for (auto model : models_) {
    model->populate_history(history);
  }
  return;
}

void CombinedInelasticity::init_history(History & history) const
{
  for (auto model : models_) {
    model->init_history(history);
  }
  return;
}

double CombinedInelasticity::strength(const History & history, Lattice & L,
                                      double T) const
{
  // May need to think on this
  double sum = 0.0;
  for (auto model : models_) {
    sum += model->strength(history, L, T);
  }
  return sum;
}

Symmetric CombinedInelasticity::d_p(const Symmetric & stress, const Orientation & Q,
                              const History & history,
                              Lattice & lattice, double T) const
{
  Symmetric sum;
  for (auto model : models_) {
    sum += model->d_p(stress, Q, history, lattice, T);
  }
  return sum;
}

SymSymR4 CombinedInelasticity::d_d_p_d_stress(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  SymSymR4 sum;
  for (auto model : models_) {
    sum += model->d_d_p_d_stress(stress, Q, history, lattice, T);
  }
  return sum;
}

History CombinedInelasticity::d_d_p_d_history(
    const Symmetric & stress, const Orientation & Q,
    const History & history,
    Lattice & lattice, double T) const
{
  History hist;
  for (auto model : models_) {
    hist.add_union(model->d_d_p_d_history(stress, Q, history, lattice, T));
  }
  return hist;
}

History CombinedInelasticity::history_rate(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history,
                                     Lattice & lattice, double T) const
{
  History hist;
  for (auto model : models_) {
    hist.add_union(model->history_rate(stress, Q, history, lattice, T));
  }
  return hist;
}

History CombinedInelasticity::d_history_rate_d_stress(const Symmetric & stress, 
                                                const Orientation & Q,
                                                const History & history,
                                                Lattice & lattice, double T) const
{
  History hist;
  for (auto model : models_) {
    hist.add_union(model->d_history_rate_d_stress(stress, Q, history, lattice, T));
  }
  return hist;
}

History CombinedInelasticity::d_history_rate_d_history(const Symmetric & stress,
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & lattice, double T) const
{
  History hist;
  for (auto model : models_) {
    hist.add_union(model->d_history_rate_d_history(stress, Q, history, lattice, T));
  }
  return hist;
}

Skew CombinedInelasticity::w_p(const Symmetric & stress, const Orientation & Q,
                         const History & history,
                         Lattice & lattice, double T) const
{
  Skew sum;
  for (auto model : models_) {
    sum += model->w_p(stress, Q, history, lattice, T);
  }
  return sum;
}

SkewSymR4 CombinedInelasticity::d_w_p_d_stress(const Symmetric & stress, 
                                       const Orientation & Q,
                                       const History & history,
                                       Lattice & lattice, double T) const
{
  SkewSymR4 sum;
  for (auto model : models_) {
    sum += model->d_w_p_d_stress(stress, Q, history, lattice, T);
  }
  return sum;
}

History CombinedInelasticity::d_w_p_d_history(const Symmetric & stress,
                                        const Orientation & Q,
                                        const History & history,
                                        Lattice & lattice, double T) const
{
  History hist;
  for (auto model : models_) {
    hist.add_union(model->d_w_p_d_history(stress, Q, history, lattice, T));
  }
  return hist;
}

} // namespace neml
