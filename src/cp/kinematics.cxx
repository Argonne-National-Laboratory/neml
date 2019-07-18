#include "kinematics.h"

namespace neml {

SymSym KinematicModel::d_stress_rate_d_d_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  return SymSym();
}

SymSkew KinematicModel::d_stress_rate_d_w_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  return SymSkew();
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

History StandardKinematicModel::decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T)
{
  History constant;

  constant.add<Skew>("espin");
  constant.add<SymSym>("C");
  constant.add<SymSym>("S");

  constant.get<Skew>("espin") = spin(stress, d, w, Q, history, lattice, T);
  constant.get<SymSym>("C") = emodel_->C(T,Q);
  constant.get<SymSym>("S") = emodel_->S(T,Q);

  return constant;
}

Symmetric StandardKinematicModel::stress_rate(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  Symmetric e = fixed.get<SymSym>("S").dot(stress);

  Skew O_s = fixed.get<Skew>("espin") + imodel_->w_p(stress, Q, history, lattice, T);

  Symmetric dp = imodel_->d_p(stress, Q, history, lattice, T);

  Symmetric net = Symmetric(e*O_s - O_s*e);

  return fixed.get<SymSym>("C").dot(d - dp - net);
}

SymSym StandardKinematicModel::d_stress_rate_d_stress(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  Symmetric e = fixed.get<SymSym>("S").dot(stress);

  Skew O_s = fixed.get<Skew>("espin") + imodel_->w_p(stress, Q, history, lattice, T);

  SymSym D1 = imodel_->d_d_p_d_stress(stress, Q, history, lattice, T);
  SymSym D2 = SymSymSkew_SkewSymSym(fixed.get<SymSym>("S"), O_s);

  SkewSym DW = imodel_->d_w_p_d_stress(stress, Q, history, lattice, T);

  SymSym D3 = SymSkewSym_SkewSymSym(DW, e);

  return -fixed.get<SymSym>("C") * (D1 + D2 + D3);
}

SymSym StandardKinematicModel::d_stress_rate_d_d(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return fixed.get<SymSym>("C");
}

SymSkew StandardKinematicModel::d_stress_rate_d_w(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return SymSkew();
}

History StandardKinematicModel::d_stress_rate_d_history(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History res = history.derivative<Symmetric>();
  History dD = imodel_->d_d_p_d_history(stress, Q, history, lattice, T);
  History dW = imodel_->d_w_p_d_history(stress, Q, history, lattice, T);

  Symmetric e = fixed.get<SymSym>("S").dot(stress);

  for (auto hvar : history.items()) {
    res.get<Symmetric>(hvar) = -fixed.get<SymSym>("C") * (
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
  return imodel_->history_rate(stress, Q, history, lattice, T);
}

History StandardKinematicModel::d_history_rate_d_stress(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  return imodel_->d_history_rate_d_stress(stress, Q, history, lattice, T);
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
  return imodel_->d_history_rate_d_history(stress, Q, history, lattice, T);
}

SymSkew StandardKinematicModel::d_stress_rate_d_w_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  Symmetric e = fixed.get<SymSym>("S").dot(stress);
  
  return -2.0 * SpecialSymSymSym(fixed.get<SymSym>("C"), e);
}

Skew StandardKinematicModel::spin(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T) const
{
  SymSym S = emodel_->S(T, Q);
  Symmetric e = S.dot(stress);

  Skew wp = imodel_->w_p(stress, Q, history, lattice, T);
  Symmetric dp = imodel_->d_p(stress, Q, history, lattice, T);
  
  Skew net = Skew(e * dp - dp * e);

  return w - wp - net;
}

Symmetric StandardKinematicModel::elastic_strains(
    const Symmetric & stress, const Orientation & Q,
    const History & history, double T)
{
  SymSym S = emodel_->S(T,Q);
  return S.dot(stress);
}

} // namespace neml
