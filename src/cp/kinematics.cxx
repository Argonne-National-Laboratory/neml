#include "kinematics.h"

namespace neml {

SymSymR4 KinematicModel::d_stress_rate_d_d_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  return SymSymR4();
}

SymSkewR4 KinematicModel::d_stress_rate_d_w_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  return SymSkewR4();
}

History KinematicModel::d_history_rate_d_d_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  return history.derivative<Symmetric>();;
}

History KinematicModel::d_history_rate_d_w_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  return history.derivative<Skew>();
}

bool KinematicModel::use_nye() const
{
  return false;
}

StandardKinematicModel::StandardKinematicModel(
    std::shared_ptr<LinearElasticModel> emodel,
    std::shared_ptr<InelasticModel> imodel) :
      emodel_(emodel), imodel_(imodel)
{

}

StandardKinematicModel::~StandardKinematicModel()
{

}

std::string StandardKinematicModel::type()
{
  return "StandardKinematicModel";
}

std::unique_ptr<NEMLObject> StandardKinematicModel::initialize(ParameterSet & params)
{
  return neml::make_unique<StandardKinematicModel>(
      params.get_object_parameter<LinearElasticModel>("emodel"),
      params.get_object_parameter<InelasticModel>("imodel"));
}

ParameterSet StandardKinematicModel::parameters()
{
  ParameterSet pset(StandardKinematicModel::type());
  
  pset.add_parameter<NEMLObject>("emodel");
  pset.add_parameter<NEMLObject>("imodel");

  return pset;
}

void StandardKinematicModel::populate_history(History & history) const
{
  imodel_->populate_history(history);
}

void StandardKinematicModel::init_history(History & history) const
{
  imodel_->init_history(history);
}

double StandardKinematicModel::strength(const History & history, Lattice & L,
                                        double T, const History & fixed) const
{
  return imodel_->strength(history, L, T, fixed);
}

History StandardKinematicModel::decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  History constant;

  constant.add<Skew>("espin");
  constant.add<SymSymR4>("C");
  constant.add<SymSymR4>("S");

  constant.get<Skew>("espin") = spin(stress, d, w, Q, history, lattice, T, fixed);
  constant.get<SymSymR4>("C") = emodel_->C(T,Q);
  constant.get<SymSymR4>("S") = emodel_->S(T,Q);
  
  // Return with whatever extra multiphysics variables you need
  return constant.add_union(fixed);
}

Symmetric StandardKinematicModel::stress_rate(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  Symmetric e = fixed.get<SymSymR4>("S").dot(stress);

  Skew O_s = fixed.get<Skew>("espin") + imodel_->w_p(stress, Q, history, lattice, T, fixed);

  Symmetric dp = imodel_->d_p(stress, Q, history, lattice, T, fixed);

  Symmetric net = Symmetric(e*O_s - O_s*e);

  return fixed.get<SymSymR4>("C").dot(d - dp - net);
}

SymSymR4 StandardKinematicModel::d_stress_rate_d_stress(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  Symmetric e = fixed.get<SymSymR4>("S").dot(stress);

  Skew O_s = fixed.get<Skew>("espin") + imodel_->w_p(stress, Q, history, lattice, T, fixed);

  SymSymR4 D1 = imodel_->d_d_p_d_stress(stress, Q, history, lattice, T, fixed);
  SymSymR4 D2 = SymSymR4Skew_SkewSymR4SymR4(fixed.get<SymSymR4>("S"), O_s);

  SkewSymR4 DW = imodel_->d_w_p_d_stress(stress, Q, history, lattice, T, fixed);

  SymSymR4 D3 = SymSkewR4Sym_SkewSymR4SymR4(DW, e);

  return -fixed.get<SymSymR4>("C") * (D1 + D2 + D3);
}

SymSymR4 StandardKinematicModel::d_stress_rate_d_d(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return fixed.get<SymSymR4>("C");
}

SymSkewR4 StandardKinematicModel::d_stress_rate_d_w(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return SymSkewR4();
}

History StandardKinematicModel::d_stress_rate_d_history(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History res = history.derivative<Symmetric>();
  History dD = imodel_->d_d_p_d_history(stress, Q, history, lattice, T, fixed);
  History dW = imodel_->d_w_p_d_history(stress, Q, history, lattice, T, fixed);

  Symmetric e = fixed.get<SymSymR4>("S").dot(stress);

  for (auto hvar : history.items()) {
    res.get<Symmetric>(hvar) = -fixed.get<SymSymR4>("C") * (
        dD.get<Symmetric>(hvar) + Symmetric(e*dW.get<Skew>(hvar) - dW.get<Skew>(hvar) * e));
  }

  return res;
}

History StandardKinematicModel::history_rate(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return imodel_->history_rate(stress, Q, history, lattice, T, fixed);
}

History StandardKinematicModel::d_history_rate_d_stress(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return imodel_->d_history_rate_d_stress(stress, Q, history, lattice, T, fixed);
}

History StandardKinematicModel::d_history_rate_d_d(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return history.derivative<Symmetric>();
}

History StandardKinematicModel::d_history_rate_d_w(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return history.derivative<Skew>();
}

History StandardKinematicModel::d_history_rate_d_history(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return imodel_->d_history_rate_d_history(stress, Q, history, lattice, T, fixed);
}

SymSkewR4 StandardKinematicModel::d_stress_rate_d_w_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  Symmetric e = fixed.get<SymSymR4>("S").dot(stress);
  
  return -2.0 * SpecialSymSymR4Sym(fixed.get<SymSymR4>("C"), e);
}

Skew StandardKinematicModel::spin(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  SymSymR4 S = emodel_->S(T, Q);
  Symmetric e = S.dot(stress);

  Skew wp = imodel_->w_p(stress, Q, history, lattice, T, fixed);
  Symmetric dp = imodel_->d_p(stress, Q, history, lattice, T, fixed);
  
  Skew net = Skew(e * dp - dp * e);

  return w - wp - net;
}

Symmetric StandardKinematicModel::elastic_strains(
    const Symmetric & stress, const Orientation & Q,
    const History & history, double T)
{
  SymSymR4 S = emodel_->S(T,Q);
  return S.dot(stress);
}

bool StandardKinematicModel::use_nye() const
{
  return imodel_->use_nye();
}

} // namespace neml
