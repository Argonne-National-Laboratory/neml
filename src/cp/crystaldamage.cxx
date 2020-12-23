#include "crystaldamage.h"

#include "../math/projections.h"

namespace neml {

CrystalDamageModel::CrystalDamageModel(std::vector<std::string> vars) :
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

void CrystalDamageModel::setup(Lattice & lattice)
{
  return;
}

void CrystalDamageModel::populate_history(History & history) const
{
  for (auto name : varnames_)
    history.add<double>(name);
}

NilDamageModel::NilDamageModel() :
    CrystalDamageModel({"whatever"})
{

}

std::string NilDamageModel::type()
{
  return "NilDamageModel";
}

std::unique_ptr<NEMLObject> NilDamageModel::initialize(ParameterSet & params)
{
  return neml::make_unique<NilDamageModel>();
}

ParameterSet NilDamageModel::parameters()
{
  ParameterSet pset(NilDamageModel::type());

  return pset;
}

void NilDamageModel::init_history(History & history) const
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
                                    const SlipRule & slip, double T) const
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
                                          const SlipRule & slip, double T) const
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
                                           double T) const
{
  History res;
  res.add<double>("whatever_whatever");
  res.zero();
  return res;
}

} // namespace neml
