#include "cp/kinematics.h"

namespace neml {

KinematicModel::KinematicModel(ParameterSet & params) :
    HistoryNEMLObject(params)
{
}

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

StandardKinematicModel::StandardKinematicModel(ParameterSet & params) :
    KinematicModel(params),
    emodel_(params.get_object_parameter<LinearElasticModel>("emodel")),
    imodel_(params.get_object_parameter<InelasticModel>("imodel"))
{

}

std::string StandardKinematicModel::type()
{
  return "StandardKinematicModel";
}

std::unique_ptr<NEMLObject> StandardKinematicModel::initialize(ParameterSet & params)
{
  return neml::make_unique<StandardKinematicModel>(params);
}

ParameterSet StandardKinematicModel::parameters()
{
  ParameterSet pset(StandardKinematicModel::type());
  
  pset.add_parameter<NEMLObject>("emodel");
  pset.add_parameter<NEMLObject>("imodel");

  return pset;
}

void StandardKinematicModel::populate_hist(History & history) const
{
  imodel_->populate_hist(history);
}

void StandardKinematicModel::init_hist(History & history) const
{
  imodel_->init_hist(history);
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
    const Symmetric & stress, Lattice & lattice, const Orientation & Q,
    const History & history, double T)
{
  SymSymR4 S = emodel_->S(T,Q);
  return S.dot(stress);
}

Symmetric StandardKinematicModel::stress_increment(
    const Symmetric & stress,
    const Symmetric & D, const Skew & W, double dt, Lattice & lattice, const Orientation & Q,
    const History & history, double T)
{
  SymSymR4 C = emodel_->C(T,Q);
  return C.dot(D) * dt;
}

bool StandardKinematicModel::use_nye() const
{
  return imodel_->use_nye();
}

DamagedStandardKinematicModel::DamagedStandardKinematicModel(
    ParameterSet & params) :
      StandardKinematicModel(params),
      dmodel_(params.get_object_parameter<CrystalDamageModel>("dmodel")),
      amodel_(params.get_object_parameter<AsaroInelasticity>("imodel"))
{
}

std::string DamagedStandardKinematicModel::type()
{
  return "DamagedStandardKinematicModel";
}

std::unique_ptr<NEMLObject> DamagedStandardKinematicModel::initialize(ParameterSet & params)
{
  return neml::make_unique<DamagedStandardKinematicModel>(params);
}

ParameterSet DamagedStandardKinematicModel::parameters()
{
  ParameterSet pset(DamagedStandardKinematicModel::type());
  
  pset.add_parameter<NEMLObject>("emodel");
  pset.add_parameter<NEMLObject>("imodel");
  pset.add_parameter<NEMLObject>("dmodel");

  return pset;
}

void DamagedStandardKinematicModel::populate_hist(History & history) const
{
  imodel_->populate_hist(history);
  dmodel_->populate_hist(history);
}

void DamagedStandardKinematicModel::init_hist(History & history) const
{
  imodel_->init_hist(history);
  dmodel_->init_hist(history);
}

Symmetric DamagedStandardKinematicModel::stress_rate(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History ihist = ihist_(history);
  History dhist = dhist_(history);

  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  SymSymR4 Pi = P.inverse();
  Symmetric stress_p = Pi.dot(stress);

  Skew O_s = fixed.get<Skew>("espin") + imodel_->w_p(stress_p, Q, ihist, lattice, T, fixed);
  Symmetric dp = imodel_->d_p(stress_p, Q, ihist, lattice, T, fixed);

  Symmetric net = Symmetric(stress_p*O_s - O_s*stress_p);

  return P.dot(fixed.get<SymSymR4>("C").dot(d - dp)) - net;
}

SymSymR4 DamagedStandardKinematicModel::d_stress_rate_d_stress(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History ihist = ihist_(history);
  History dhist = dhist_(history);

  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  SymSymSymR6 Pp = dmodel_->d_projection_d_stress(stress, dhist, Q, lattice,
                                                  amodel_->slip_rule(), T);

  SymSymR4 Pi = P.inverse();
  Symmetric stress_p = Pi.dot(stress);

  Symmetric dp = imodel_->d_p(stress_p, Q, ihist, lattice, T, fixed);

  Skew O_s = fixed.get<Skew>("espin") + imodel_->w_p(stress_p, Q, ihist, lattice, T, fixed);

  SymSymR4 D1 = imodel_->d_d_p_d_stress(stress_p, Q, ihist, lattice, T,
                                        fixed);
  SymSymR4 D2 = SymSymR4Skew_SkewSymR4SymR4(SymSymR4::id(), O_s);
  
  SkewSymR4 DW = imodel_->d_w_p_d_stress(stress_p, Q, ihist, lattice, T,
                                         fixed);

  SymSymR4 D3 = SymSkewR4Sym_SkewSymR4SymR4(DW, stress_p);

  return (Pp.dot_k(d-dp) - (P.dot(fixed.get<SymSymR4>("C")) * D1 + D2 +
                            D3)).dot(Pi);
}

SymSymR4 DamagedStandardKinematicModel::d_stress_rate_d_d(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History dhist = dhist_(history);

  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  return P.dot(fixed.get<SymSymR4>("C"));
}

History DamagedStandardKinematicModel::d_stress_rate_d_history(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  // Base names
  auto inames = inames_();

  History ihist = ihist_(history);
  History dhist = dhist_(history);

  // Gets all the variables
  History res = history.derivative<Symmetric>();

  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  History pdhist = dmodel_->d_projection_d_history(stress, dhist,
                                                  Q, lattice,
                                                  amodel_->slip_rule(), T);

  SymSymR4 Pi = P.inverse();
  Symmetric stress_p = Pi.dot(stress);

  // imodel internal variables
  History dD = imodel_->d_d_p_d_history(stress_p, Q, ihist, lattice, T, fixed);
  History dW = imodel_->d_w_p_d_history(stress_p, Q, ihist, lattice, T, fixed);

  for (auto hvar : inames) {
    res.get<Symmetric>(hvar) = -P.dot(fixed.get<SymSymR4>("C")) * (
        dD.get<Symmetric>(hvar)) - Symmetric(stress_p*dW.get<Skew>(hvar) -
                                             dW.get<Skew>(hvar) * stress_p);
  }

  // dmodel internal variables
  Symmetric dp = imodel_->d_p(stress_p, Q, ihist, lattice, T, fixed);
  SymSymR4 dS = d_stress_partial(stress, d, w, Q, history, lattice, T,
                                       fixed);
  Symmetric partial = fixed.get<SymSymR4>("C").dot(d-dp);

  for (auto hvar : dhist.items()) {
    res.get<Symmetric>(hvar) =
        pdhist.get<SymSymR4>(hvar).dot(partial) -
        dS.dot(Pi.dot(pdhist.get<SymSymR4>(hvar).dot(stress_p)));
  }

  return res;
}

History DamagedStandardKinematicModel::history_rate(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History ihist = ihist_(history);
  History dhist = dhist_(history);

  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  SymSymR4 Pi = P.inverse();
  Symmetric stress_p = Pi.dot(stress);

  History base = imodel_->history_rate(stress_p, Q, ihist, lattice, T, fixed);
  base.add_union(dmodel_->damage_rate(stress_p, history, Q, lattice, 
                                      amodel_->slip_rule(), T, fixed));
  return base;
}

History DamagedStandardKinematicModel::d_history_rate_d_stress(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History ihist = ihist_(history);
  History dhist = dhist_(history);

  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  SymSymR4 Pi = P.inverse();
  Symmetric stress_p = Pi.dot(stress);

  History base = imodel_->d_history_rate_d_stress(stress_p, Q, ihist, lattice, T, fixed);
  base.add_union(dmodel_->d_damage_d_stress(stress_p, history, Q, lattice,
                                            amodel_->slip_rule(), T, fixed));
  // Need to postmultiply by Pi ...
  return base.postmultiply(Pi);
}

History DamagedStandardKinematicModel::d_history_rate_d_history(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  // We need the names of the internal variables for both the base and damage
  // models
  History ihist = ihist_(history);
  History dhist = dhist_(history);
  std::vector<std::string> tnames(ihist.items());
  tnames.insert(tnames.end(), dhist.items().begin(), dhist.items().end());

  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  SymSymR4 Pi = P.inverse();
  Symmetric stress_p = Pi.dot(stress);
  History dh = dmodel_->d_projection_d_history(stress, dhist, Q, lattice,
                                               amodel_->slip_rule(), T);

  // This assumes the *inelastic* variables are disjoint from the *damage*
  // variables but NOT the other way around
  History base =  imodel_->d_history_rate_d_history(stress_p, Q, ihist, lattice, T, fixed);
  base.add_union(dmodel_->d_damage_d_history(stress_p, history, Q, lattice,
                                            amodel_->slip_rule(), T, fixed));

  // Add zeros for the cross terms
  History id = ihist.history_derivative(dhist);
  id.zero();
  base.add_union(id);

  // There is a chain rule term...
  History partial = imodel_->d_history_rate_d_stress(stress_p, Q, ihist, lattice, T, fixed);
  partial.add_union(dmodel_->d_damage_d_stress(stress_p, history, Q, lattice,
                                            amodel_->slip_rule(), T, fixed));
  for (auto name : tnames) {
    for (auto dname : dhist.items()) {
      SymSymR4 T = -Pi.dot(dh.get<SymSymR4>(dname).dot(Pi));
      base.get<double>(name+"_"+dname) -=
          partial.get<Symmetric>(name).contract(Pi.dot(dh.get<SymSymR4>(dname).dot(stress_p)));
    }
  }

  // Out of order now...
  std::vector<std::string> names;
  for (auto iname : tnames)
  {
    for (auto dname : tnames)
    {
      names.push_back(iname+"_"+dname);
    }
  }
 
  base.reorder(names);

  return base;
}

SymSkewR4 DamagedStandardKinematicModel::d_stress_rate_d_w_decouple(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed)
{
  History dhist = dhist_(history);
  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  SymSymR4 Pi = P.inverse();
  Symmetric stress_p = Pi.dot(stress);

  return -2.0 * SpecialSymSymR4Sym(SymSymR4::id(), stress_p);
}

// Not done
Skew DamagedStandardKinematicModel::spin(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History ihist = ihist_(history);
  History dhist = dhist_(history);

  SymSymR4 S = emodel_->S(T, Q);
  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice,
                                   amodel_->slip_rule(), T);
  Symmetric e = S.dot(P.inverse()).dot(stress);

  Skew wp = imodel_->w_p(stress, Q, ihist, lattice, T, fixed);
  Symmetric dp = imodel_->d_p(stress, Q, ihist, lattice, T, fixed);
  
  Skew net = Skew(e * dp - dp * e);

  return w - wp - net;
}

Symmetric DamagedStandardKinematicModel::elastic_strains(
    const Symmetric & stress, Lattice & lattice, const Orientation & Q,
    const History & history, double T)
{
  SymSymR4 S = emodel_->S(T,Q);
  SymSymR4 P = dmodel_->projection(stress, history, Q, lattice,
                                   amodel_->slip_rule(), T);
  return S.dot(P.inverse().dot(stress));
}

Symmetric DamagedStandardKinematicModel::stress_increment(
    const Symmetric & stress,
    const Symmetric & D, const Skew & W, double dt, Lattice & lattice, const Orientation & Q,
    const History & history, double T)
{
  SymSymR4 C = emodel_->C(T,Q);
  SymSymR4 P = dmodel_->projection(stress, history, Q, lattice,
                                   amodel_->slip_rule(), T);
  return P.dot(C.dot(D*dt));
}

std::vector<std::string> DamagedStandardKinematicModel::inames_() const
{
  History ihist;
  imodel_->populate_hist(ihist);
  return ihist.items();
}

std::vector<std::string> DamagedStandardKinematicModel::dnames_() const
{
  return dmodel_->varnames();
}

History DamagedStandardKinematicModel::ihist_(const History & hist) const
{
  return hist.split(inames_(), false);
}

History DamagedStandardKinematicModel::dhist_(const History & hist) const
{
  return hist.split(inames_(), true);
}

SymSymR4 DamagedStandardKinematicModel::d_stress_partial(
    const Symmetric & stress, const Symmetric & d,
    const Skew & w, const Orientation & Q,
    const History & history, Lattice & lattice,
    double T, const History & fixed) const
{
  History ihist = ihist_(history);
  History dhist = dhist_(history);

  SymSymR4 P = dmodel_->projection(stress, dhist, Q, lattice, 
                                   amodel_->slip_rule(), T);
  SymSymR4 Pi = P.inverse();
  Symmetric stress_p = Pi.dot(stress);

  Symmetric dp = imodel_->d_p(stress_p, Q, ihist, lattice, T, fixed);

  Skew O_s = fixed.get<Skew>("espin") + imodel_->w_p(stress_p, Q, ihist, lattice, T, fixed);

  SymSymR4 D1 = imodel_->d_d_p_d_stress(stress_p, Q, ihist, lattice, T,
                                        fixed);
  SymSymR4 D2 = SymSymR4Skew_SkewSymR4SymR4(SymSymR4::id(), O_s);
  
  SkewSymR4 DW = imodel_->d_w_p_d_stress(stress_p, Q, ihist, lattice, T,
                                         fixed);

  SymSymR4 D3 = SymSkewR4Sym_SkewSymR4SymR4(DW, stress_p);

  return -(P.dot(fixed.get<SymSymR4>("C")) * D1 + D2 + D3);
}

} // namespace neml
