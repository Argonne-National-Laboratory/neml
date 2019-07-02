#include "kinematics.h"

namespace neml {

StandardKinematicModel::StandardKinematicModel(
    std::shared_ptr<LinearElasticModel> emodel,
    std::shared_ptr<InelasticModel> imodel) :
      emodel_(emodel), imodel_(imodel), espin_(Skew())
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

void StandardKinematicModel::decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T)
{
  // Lag the rotational update
  espin_ = spin(stress, d, w, Q, history, lattice, T);
}

Symmetric StandardKinematicModel::stress_rate(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{
  SymSym C = emodel_->C(T, Q);
  SymSym S = emodel_->S(T, Q);

  Symmetric e = S.dot(stress);

  Skew O_s = espin_ + imodel_->w_p(stress, Q, history, lattice, T);

  Symmetric dp = imodel_->d_p(stress, Q, history, lattice, T);

  Symmetric net = Symmetric(e*O_s - O_s*e);

  return C.dot(d - dp - net);
}

SymSym StandardKinematicModel::d_stress_rate_d_stress(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{

}

SymSym StandardKinematicModel::d_stress_rate_d_d(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{

}

SymSkew StandardKinematicModel::d_stress_rate_d_w(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{

}

History StandardKinematicModel::d_stress_rate_d_history(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{

}

History StandardKinematicModel::history_rate(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{
  return imodel_->history_rate(stress, Q, history, lattice, T);
}

History StandardKinematicModel::d_history_rate_d_stress(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{

}

History StandardKinematicModel::d_history_rate_d_d(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{

}

History StandardKinematicModel::d_history_rate_d_w(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{

}

History StandardKinematicModel::d_history_rate_d_history(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{

}

Skew StandardKinematicModel::spin(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, const Lattice & lattice,
    double T) const
{
  SymSym S = emodel_->S(T, Q);
  Symmetric e = S.dot(stress);

  Skew wp = imodel_->w_p(stress, Q, history, lattice, T);
  Symmetric dp = imodel_->d_p(stress, Q, history, lattice, T);
  
  Skew net = Skew(e * dp - dp * e);

  return w - wp - net;
}

} // namespace neml
