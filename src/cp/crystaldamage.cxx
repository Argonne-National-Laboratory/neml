#include "cp/crystaldamage.h"

#include "math/projections.h"

namespace neml {

CrystalDamageModel::CrystalDamageModel(ParameterSet & params, 
                                       std::vector<std::string> vars) :
    HistoryNEMLObject(params),
    varnames_(vars)
{

}

size_t CrystalDamageModel::nvars() const
{
  return varnames_.size();
}

std::vector<std::string> CrystalDamageModel::varnames() const
{
  return varnames_;
}

void CrystalDamageModel::set_varnames(std::vector<std::string> names)
{
  varnames_ = names;
}

void CrystalDamageModel::populate_hist(History & history) const
{
  for (auto name : varnames_)
    history.add<double>(name);
}

NilDamageModel::NilDamageModel(ParameterSet & params) :
    CrystalDamageModel(params, {"whatever"})
{

}

std::string NilDamageModel::type()
{
  return "NilDamageModel";
}

std::unique_ptr<NEMLObject> NilDamageModel::initialize(ParameterSet & params)
{
  return neml::make_unique<NilDamageModel>(params);
}

ParameterSet NilDamageModel::parameters()
{
  ParameterSet pset(NilDamageModel::type());

  return pset;
}

void NilDamageModel::init_hist(History & history) const
{
  history.get<double>("whatever") = 0.0;
}

SymSymR4 NilDamageModel::projection(const Symmetric & stress,
                                    const History & damage, 
                                    const Orientation & Q, Lattice & lattice,
                                    const SlipRule & slip, double T)
{
  return SymSymR4::id();
}

SymSymSymR6 NilDamageModel::d_projection_d_stress(const Symmetric & stress,
                                                  const History & damage,
                                                  const Orientation & Q,
                                                  Lattice & lattice, 
                                                  const SlipRule & slip,
                                                  double T)
{
  SymSymSymR6 zero;
  return zero;
}

History NilDamageModel::d_projection_d_history(const Symmetric & stress,
                                               const History & damage, 
                                               const Orientation & Q,
                                               Lattice & lattice, 
                                               const SlipRule & slip, double T)
{
  History res;
  res.add<SymSymR4>("whatever");
  res.zero();
  return res;
}

History NilDamageModel::damage_rate(const Symmetric & stress,
                                    const History & history, 
                                    const Orientation & Q, Lattice & lattice, 
                                    const SlipRule & slip, double T,
                                    const History & fixed) const
{
  History res;
  res.add<double>("whatever");
  res.get<double>("whatever") = 0;
  return res;
}

/// Derivative of each damage with respect to stress
History NilDamageModel::d_damage_d_stress(const Symmetric & stress,
                                          const History & history, 
                                          const Orientation & Q,
                                          Lattice & lattice, 
                                          const SlipRule & slip, double T,
                                          const History & fixed) const
{
  History res;
  res.add<Symmetric>("whatever");
  res.zero();
  return res;
}

/// Derivative of damage with respect to history
History NilDamageModel::d_damage_d_history(const Symmetric & stress,
                                           const History & history,
                                           const Orientation & Q,
                                           Lattice & lattice,
                                           const SlipRule & slip, 
                                           double T, 
                                           const History & fixed) const
{
  History dvars = history.subset(varnames_);
  History res = dvars.history_derivative(history);
  
  res.zero();

  return res;
}

PlanarDamageModel::PlanarDamageModel(ParameterSet & params) :
      CrystalDamageModel(params, {}), 
      damage_(params.get_object_parameter<SlipPlaneDamage>("damage")), 
      shear_transform_(params.get_object_parameter<TransformationFunction>("shear_transform")),
      normal_transform_(params.get_object_parameter<TransformationFunction>("normal_transform")), 
      lattice_(params.get_object_parameter<Lattice>("lattice"))
{
  std::vector<std::string> varnames;
  for (size_t i = 0; i < lattice_->nplanes(); i++)
    varnames.push_back("slip_damage_"+std::to_string(i));
  set_varnames(varnames);
}

std::string PlanarDamageModel::type()
{
  return "PlanarDamageModel";
}

std::unique_ptr<NEMLObject> PlanarDamageModel::initialize(ParameterSet & params)
{
  return neml::make_unique<PlanarDamageModel>(params);
}

ParameterSet PlanarDamageModel::parameters()
{
  ParameterSet pset(PlanarDamageModel::type());

  pset.add_parameter<NEMLObject>("damage");
  pset.add_parameter<NEMLObject>("shear_transform");
  pset.add_parameter<NEMLObject>("normal_transform");
  pset.add_parameter<NEMLObject>("lattice");

  return pset;
}

void PlanarDamageModel::init_hist(History & history) const
{
  for (auto vn : varnames_)
    history.get<double>(vn) = damage_->setup();
}

SymSymR4 PlanarDamageModel::projection(const Symmetric & stress,
                                    const History & damage, 
                                    const Orientation & Q, Lattice & lattice,
                                    const SlipRule & slip, double T)
{
  SymSymR4 P = SymSymR4::id();
  for (size_t i = 0; i < lattice.nplanes(); i++) {
    Vector n = Q.apply(lattice.unique_planes()[i]); // ...
    
    SymSymR4 P_s = shear_projection_ss(n);
    SymSymR4 P_n = normal_projection_ss(n);
    
    double ns = stress.dot(n).dot(n);
    double d = damage.get<double>(varnames_[i]);
    double ds = shear_transform_->map(d, ns);
    double dn = normal_transform_->map(d, ns);

    P = P.dot(SymSymR4::id() - ds * P_s - dn * P_n);
  }

  return P;
}

SymSymSymR6 PlanarDamageModel::d_projection_d_stress(const Symmetric & stress,
                                                  const History & damage,
                                                  const Orientation & Q,
                                                  Lattice & lattice, 
                                                  const SlipRule & slip,
                                                  double T)
{
  SymSymSymR6 Pd;
  
  for (size_t i = 0; i < lattice.nplanes(); i++) {
    SymSymR4 before = SymSymR4::id();
    SymSymSymR6 middle;
    SymSymR4 after = SymSymR4::id();
    for (size_t j = 0; j < lattice.nplanes(); j++) {
      Vector n = Q.apply(lattice.unique_planes()[j]); // ...
      
      SymSymR4 P_s = shear_projection_ss(n);
      SymSymR4 P_n = normal_projection_ss(n);
      
      double ns = stress.dot(n).dot(n);
      double d = damage.get<double>(varnames_[j]);
      double ds = shear_transform_->map(d, ns);
      double dn = normal_transform_->map(d, ns);
      
      if (i < j) {
        before = before.dot(SymSymR4::id() - ds * P_s - dn * P_n);
      }
      else if (i == j) {
        middle = -outer_product_k(P_s, 
                                  shear_transform_->d_map_d_normal(d, ns) *
                                  Symmetric(n.outer(n)))
            - outer_product_k(P_n,
                              normal_transform_->d_map_d_normal(d, ns) *
                              Symmetric(n.outer(n)));
      }
      else {
        after = after.dot(SymSymR4::id() - ds * P_s - dn * P_n);
      }
    }
    Pd += middle.middle_dot_after(after).middle_dot_before(before);
  }

  return Pd;
}

History PlanarDamageModel::d_projection_d_history(const Symmetric & stress,
                                               const History & damage, 
                                               const Orientation & Q,
                                               Lattice & lattice, 
                                               const SlipRule & slip, double T)
{
  History res;
  
  for (size_t i = 0; i < lattice.nplanes(); i++) {
    res.add<SymSymR4>(varnames_[i]);
    res.get<SymSymR4>(varnames_[i]) = SymSymR4::id();
    for (size_t j = 0; j < lattice.nplanes(); j++) {
      Vector n = Q.apply(lattice.unique_planes()[j]); // ...
      
      SymSymR4 P_s = shear_projection_ss(n);
      SymSymR4 P_n = normal_projection_ss(n);
      
      double ns = stress.dot(n).dot(n);
      double d = damage.get<double>(varnames_[j]);
      double ds = shear_transform_->map(d, ns);
      double dn = normal_transform_->map(d, ns);
      
      if (i != j) {
        res.get<SymSymR4>(varnames_[i]) = 
            res.get<SymSymR4>(varnames_[i]).dot(SymSymR4::id() - 
                                                ds * P_s - dn * P_n);
      }
      else {
        res.get<SymSymR4>(varnames_[i]) = res.get<SymSymR4>(varnames_[i]).dot(
            -P_s * shear_transform_->d_map_d_damage(d, ns) 
            - P_n * normal_transform_->d_map_d_damage(d, ns));
      }
    }
  }

  return res;
}

History PlanarDamageModel::damage_rate(
    const Symmetric & stress, const History & history, 
    const Orientation & Q, Lattice & lattice, const SlipRule & slip,
    double T, const History & fixed) const
{
  History imodel = inelastic_history_(history);

  History res;

  for (size_t i = 0; i < lattice.nplanes(); i++) {
    res.add<double>(varnames_[i]);

    Vector n = Q.apply(lattice.unique_planes()[i]);
    auto systems = lattice.plane_systems(i);
    std::vector<double> taus(systems.size());
    std::vector<double> gammas(systems.size());

    for (size_t k = 0; k < systems.size(); k++) {
      size_t g = systems[k].first;
      size_t ii = systems[k].second;

      taus[k] = lattice.shear(g, ii, Q, stress);
      gammas[k] = slip.slip(g, ii, stress, Q, imodel, lattice, T, fixed);
    }
    
    double ns = stress.dot(n).dot(n);
    double d = history.get<double>(varnames_[i]);

    res.get<double>(varnames_[i]) = damage_->damage_rate(taus, gammas, ns, d);
  }
  return res;
}

/// Derivative of each damage with respect to stress
History PlanarDamageModel::d_damage_d_stress(const Symmetric & stress,
                                          const History & history, 
                                          const Orientation & Q,
                                          Lattice & lattice, 
                                          const SlipRule & slip, double T,
                                          const History & fixed) const
{
  History imodel = inelastic_history_(history);

  History res = history.subset(varnames_).derivative<Symmetric>();
  
  for (size_t i = 0; i < lattice.nplanes(); i++) {
    Vector n = Q.apply(lattice.unique_planes()[i]);
    auto systems = lattice.plane_systems(i);
    std::vector<double> taus(systems.size());
    std::vector<double> gammas(systems.size());

    for (size_t k = 0; k < systems.size(); k++) {
      size_t g = systems[k].first;
      size_t ii = systems[k].second;

      taus[k] = lattice.shear(g, ii, Q, stress);
      gammas[k] = slip.slip(g, ii, stress, Q, imodel, lattice, T, fixed);
    }
    
    double ns = stress.dot(n).dot(n);
    double d = history.get<double>(varnames_[i]);

    Symmetric ndir = Symmetric(n.outer(n));

    res.get<Symmetric>(varnames_[i]) = 
        damage_->d_damage_rate_d_normal(taus, gammas, ns, d) * ndir;
    
    auto dshear = damage_->d_damage_rate_d_shear(taus, gammas, ns, d);
    auto dslip  = damage_->d_damage_rate_d_slip(taus, gammas, ns, d);

    for (size_t k = 0; k < systems.size(); k++) {
      size_t g = systems[k].first;
      size_t ii = systems[k].second;
      res.get<Symmetric>(varnames_[i]) += 
          dshear[k] * lattice.d_shear(g, ii, Q, stress)
        + dslip[k] * slip.d_slip_d_s(g, ii, stress, Q, imodel, lattice,
                                     T, fixed);
    }
  }

  return res;
}

/// Derivative of damage with respect to history
History PlanarDamageModel::d_damage_d_history(const Symmetric & stress,
                                           const History & history,
                                           const Orientation & Q,
                                           Lattice & lattice,
                                           const SlipRule & slip, 
                                           double T,
                                           const History & fixed) const
{
  History imodel = inelastic_history_(history);
  History dvars = damage_history_(history);
  History res = dvars.history_derivative(history);
  res.zero();

  for (size_t i = 0; i < lattice.nplanes(); i++) {
    Vector n = Q.apply(lattice.unique_planes()[i]);
    auto systems = lattice.plane_systems(i);
    std::vector<double> taus(systems.size());
    std::vector<double> gammas(systems.size());

    for (size_t k = 0; k < systems.size(); k++) {
      size_t g = systems[k].first;
      size_t ii = systems[k].second;

      taus[k] = lattice.shear(g, ii, Q, stress);
      gammas[k] = slip.slip(g, ii, stress, Q, imodel, lattice, T, fixed);
    }
    
    double ns = stress.dot(n).dot(n);
    double d = history.get<double>(varnames_[i]);

    Symmetric ndir = Symmetric(n.outer(n));

    // The damage model terms arise only from the object itself
    res.get<double>(varnames_[i]+"_"+varnames_[i]) = 
        damage_->d_damage_rate_d_damage(taus, gammas, ns, d);
    
    // The cross terms come through the gamma_dot dependence
    auto dg = damage_->d_damage_rate_d_slip(taus, gammas, ns, d);

    for (size_t k = 0; k < systems.size(); k++) {
      size_t g = systems[k].first;
      size_t ii = systems[k].second;
      History sd = slip.d_slip_d_h(g, ii, stress, Q, imodel, lattice, 
                                   T, fixed);
      for (auto other : sd.items()) {
        res.get<double>(varnames_[i]+"_"+other) += dg[k] *
            sd.get<double>(other);
      }
    }
  }

  return res;
}

History PlanarDamageModel::damage_history_(const History & total) const
{
  return total.subset(varnames_);
}

History PlanarDamageModel::inelastic_history_(const History & total) const
{
  auto names = total.items();
  for (auto name : varnames_) {
    names.erase(std::remove(names.begin(), names.end(), name), names.end());
  }

  return total.subset(names);
}

SlipPlaneDamage::SlipPlaneDamage(ParameterSet & params) :
    NEMLObject(params)
{

}

WorkPlaneDamage::WorkPlaneDamage(ParameterSet & params) :
    SlipPlaneDamage(params)
{

}

std::string WorkPlaneDamage::type()
{
  return "WorkPlaneDamage";
}

std::unique_ptr<NEMLObject> WorkPlaneDamage::initialize(ParameterSet & params)
{
  return neml::make_unique<WorkPlaneDamage>(params);
}

ParameterSet WorkPlaneDamage::parameters()
{
  ParameterSet pset(WorkPlaneDamage::type());

  return pset;
}

double WorkPlaneDamage::setup() const
{
  return 0;
}

double WorkPlaneDamage::damage_rate(
    const std::vector<double> & shears, const std::vector<double> & sliprates,
    double normal_stress, double damage)
{
  double rate = 0;
  for (size_t i = 0; i < shears.size(); i++) 
    rate += shears[i] * sliprates[i];

  return rate;
}

std::vector<double> WorkPlaneDamage::d_damage_rate_d_shear(
    const std::vector<double> & shears, const std::vector<double> & sliprates,
    double normal_stress, double damage)
{
  std::vector<double> deriv;
  deriv.resize(shears.size());
  for (size_t i = 0; i < shears.size(); i++)
    deriv[i] = sliprates[i];
  return deriv;
}

std::vector<double> WorkPlaneDamage::d_damage_rate_d_slip(
    const std::vector<double> & shears, const std::vector<double> & sliprates,
    double normal_stress, double damage)
{
  std::vector<double> deriv;
  deriv.resize(shears.size());
  for (size_t i = 0; i < shears.size(); i++)
    deriv[i] = shears[i];
  return deriv;
}

double WorkPlaneDamage::d_damage_rate_d_normal(
    const std::vector<double> & shears, const std::vector<double> & sliprates,
    double normal_stress, double damage)
{
  return 0;
}

double WorkPlaneDamage::d_damage_rate_d_damage(
    const std::vector<double> & shears, const std::vector<double> & sliprates,
    double normal_stress, double damage)
{
  return 0;
}

TransformationFunction::TransformationFunction(ParameterSet & params) :
    NEMLObject(params)
{

}

SigmoidTransformation::SigmoidTransformation(ParameterSet & params) :
    TransformationFunction(params),
    c_(params.get_parameter<double>("c")), 
    beta_(params.get_parameter<double>("beta")), 
    cut_(params.get_parameter<double>("cut"))
{

}

std::string SigmoidTransformation::type()
{
  return "SigmoidTransformation";
}

std::unique_ptr<NEMLObject> SigmoidTransformation::initialize(
    ParameterSet & params)
{
  return neml::make_unique<SigmoidTransformation>(params);
}

ParameterSet SigmoidTransformation::parameters()
{
  ParameterSet pset(SigmoidTransformation::type());
  
  pset.add_parameter<double>("c");
  pset.add_parameter<double>("beta");
  pset.add_optional_parameter<double>("cut", 1.0);

  return pset;
}

double SigmoidTransformation::map(double damage, double normal_stress)
{
  if (damage < 0)
    return 0;
  else if (damage < c_)
  {
    double dtrial = 1.0/(1.0 + std::pow(c_ / damage - 1.0, beta_));
    if (dtrial > cut_) 
      return cut_;
    else
      return dtrial;
  }
  return cut_;
}

double SigmoidTransformation::d_map_d_damage(double damage, double normal_stress)
{
  if (damage < 0)
    return 0;
  else if (damage < c_) {
    double dtrial = 1.0/(1.0 + std::pow(c_ / damage - 1.0, beta_));
    if (dtrial > cut_)
      return 0;
    else
      return beta_*c_*std::pow(damage,beta_-1)*std::pow(1.0/(c_-damage),beta_+1)
          / pow(std::pow(damage/(c_-damage),beta_)+1.0,2.0);
  }
  return 0;
}

double SigmoidTransformation::d_map_d_normal(double damage, double normal_stress)
{
  return 0;
}

SwitchTransformation::SwitchTransformation(ParameterSet & params) :
    TransformationFunction(params),
    base_(params.get_object_parameter<TransformationFunction>("base"))
{

}

std::string SwitchTransformation::type()
{
  return "SwitchTransformation";
}

std::unique_ptr<NEMLObject> SwitchTransformation::initialize(
    ParameterSet & params)
{
  return neml::make_unique<SwitchTransformation>(params);
}

ParameterSet SwitchTransformation::parameters()
{
  ParameterSet pset(SwitchTransformation::type());
  
  pset.add_parameter<NEMLObject>("base");

  return pset;
}

double SwitchTransformation::map(double damage, double normal_stress)
{
  if (normal_stress >= 0.0)
    return base_->map(damage, normal_stress);
  return 0;
}

double SwitchTransformation::d_map_d_damage(double damage, double normal_stress)
{
  if (normal_stress >= 0.0)
    return base_->d_map_d_damage(damage, normal_stress);
  return 0;
}

double SwitchTransformation::d_map_d_normal(double damage, double normal_stress)
{
  return 0;
}


} // namespace neml
