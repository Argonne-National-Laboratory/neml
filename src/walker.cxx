#include "walker.h"

namespace neml {

WalkerKremplSwitchRule::WalkerKremplSwitchRule(ParameterSet & params) :
    GeneralFlowRule(params),
    elastic_(params.get_object_parameter<LinearElasticModel>("elastic")),
    flow_(params.get_object_parameter<ViscoPlasticFlowRule>("flow")),
    lambda_(params.get_object_parameter<Interpolate>("lambda")),
    eps0_(params.get_parameter<double>("eps_ref"))
{
  cache_history_();
}

std::string WalkerKremplSwitchRule::type()
{
  return "WalkerKremplSwitchRule";
}

ParameterSet WalkerKremplSwitchRule::parameters()
{
  ParameterSet pset(WalkerKremplSwitchRule::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("flow");
  pset.add_parameter<NEMLObject>("lambda");
  pset.add_parameter<double>("eps_ref");

  return pset;
}

std::unique_ptr<NEMLObject> WalkerKremplSwitchRule::initialize(ParameterSet & params)
{
  return neml::make_unique<WalkerKremplSwitchRule>(params); 
}


void WalkerKremplSwitchRule::populate_hist(History & hist) const
{
  flow_->set_variable_prefix(get_variable_prefix());
  return flow_->populate_hist(hist);
}

void WalkerKremplSwitchRule::init_hist(History & hist) const
{
  return flow_->init_hist(hist);
}

Symmetric WalkerKremplSwitchRule::s(const Symmetric & s, const History & alpha, 
                                    const Symmetric & edot, double T, double Tdot)
{
  Symmetric erate = edot;

  Symmetric dir;
  double yv;
  flow_->g(s.data(), alpha.rawptr(), T, dir.s());
  flow_->y(s.data(), alpha.rawptr(), T, yv);
  
  double kap = kappa(edot, T);

  erate -= yv * kap * dir;

  return elastic_->C(T).dot(erate);
}

SymSymR4 WalkerKremplSwitchRule::ds_ds(const Symmetric & s, const History & alpha, 
                                       const Symmetric & edot, double T, double Tdot)
{
  double yv;
  flow_->y(s.data(), alpha.rawptr(), T, yv);

  double kap = kappa(edot, T);

  SymSymR4 work;
  flow_->dg_ds(s.data(), alpha.rawptr(), T, work.s());
  work *= -yv * kap;

  Symmetric t1;
  flow_->g(s.data(), alpha.rawptr(), T, t1.s());
  Symmetric t2;
  flow_->dy_ds(s.data(), alpha.rawptr(), T, t2.s());
  t2 *= kap;
  work -= douter(t1, t2);
  
  return elastic_->C(T).dot(work);
}

History WalkerKremplSwitchRule::ds_da(const Symmetric & s, const History & alpha, 
                                      const Symmetric & edot, double T, double Tdot)
{
  double yv;
  flow_->y(s.data(), alpha.rawptr(), T, yv);
 
  double kap = kappa(edot, T);

  History work = gather_blank_history_().derivative<Symmetric>();
  flow_->dg_da(s.data(), alpha.rawptr(), T, work.rawptr());
  work *= -yv * kap;

  Symmetric t1;
  flow_->g(s.data(), alpha.rawptr(), T, t1.s());
  History t2 = gather_blank_history_();
  flow_->dy_da(s.data(), alpha.rawptr(), T, t2.rawptr());
  t2 *= kap;
  outer_update_minus(t1.data(), 6, t2.rawptr(), nhist(), work.rawptr());
  
  return work.premultiply(elastic_->C(T));
}

SymSymR4 WalkerKremplSwitchRule::ds_de(const Symmetric & s, const History & alpha, 
                                       const Symmetric & edot, double T, double Tdot)
{
  double yv;
  flow_->y(s.data(), alpha.rawptr(), T, yv);
  
  Symmetric dkap = dkappa(edot, T);

  Symmetric g;
  flow_->g(s.data(), alpha.rawptr(), T, g.s());
  g *= yv;
  
  return elastic_->C(T).dot(SymSymR4::id() - douter(g, dkap));
}

History WalkerKremplSwitchRule::a(const Symmetric & s, const History & alpha, 
                                  const Symmetric & edot, double T, double Tdot)
{
  double dg;
  flow_->y(s.data(), alpha.rawptr(), T, dg);

  double kap = kappa(edot, T);
  
  History adot = gather_blank_history_();
  flow_->h(s.data(), alpha.rawptr(), T, adot.rawptr());
  adot *= dg * kap;
  
  History temp = gather_blank_history_();
  flow_->h_temp(s.data(), alpha.rawptr(), T, temp.rawptr());
  adot += temp * Tdot;
  flow_->h_time(s.data(), alpha.rawptr(), T, temp.rawptr());

  return adot += temp * kap;
}

History WalkerKremplSwitchRule::da_ds(const Symmetric & s, const History & alpha, 
                                      const Symmetric & edot, double T, double Tdot)
{
  double dg;
  flow_->y(s.data(), alpha.rawptr(), T, dg);

  double kap = kappa(edot, T);

  History d_adot = gather_blank_history_().derivative<Symmetric>();
  flow_->dh_ds(s.data(), alpha.rawptr(), T, d_adot.rawptr());
  d_adot *= dg * kap;
  
  History t1 = gather_blank_history_();
  flow_->h(s.data(), alpha.rawptr(), T, t1.rawptr());

  Symmetric t2;
  flow_->dy_ds(s.data(), alpha.rawptr(), T, t2.s());
  t2 *= kap;

  outer_update(t1.rawptr(), nhist(), t2.data(), 6, d_adot.rawptr());
  
  History t3 = gather_blank_history_().derivative<Symmetric>();
  flow_->dh_ds_temp(s.data(), alpha.rawptr(), T, t3.rawptr());
  d_adot += t3 * Tdot;
  flow_->dh_ds_time(s.data(), alpha.rawptr(), T, t3.rawptr());
  return d_adot + t3 * kap;
}

History WalkerKremplSwitchRule::da_da(const Symmetric & s, const History & alpha, 
                                      const Symmetric & edot, double T, double Tdot)
{
  double dg;
  flow_->y(s.data(), alpha.rawptr(), T, dg);

  double kap = kappa(edot, T);

  int nh = nhist();
  
  History d_adot = gather_blank_history_().derivative<History>();
  flow_->dh_da(s.data(), alpha.rawptr(), T, d_adot.rawptr());
  d_adot *= dg * kap;
 
  History t1 = gather_blank_history_();
  flow_->h(s.data(), alpha.rawptr(), T, t1.rawptr());
  
  History t2 = gather_blank_history_();
  flow_->dy_da(s.data(), alpha.rawptr(), T, t2.rawptr());
  t2 *= kap;

  outer_update(t1.rawptr(), nh, t2.rawptr(), nh, d_adot.rawptr());
  
  History t3 = gather_blank_history_().derivative<History>();
  flow_->dh_da_temp(s.data(), alpha.rawptr(), T, t3.rawptr());
  d_adot += t3 * Tdot;
  flow_->dh_da_time(s.data(), alpha.rawptr(), T, t3.rawptr());
  return d_adot + t3 * kap;
}

History WalkerKremplSwitchRule::da_de(const Symmetric & s, const History & alpha, 
                                      const Symmetric & edot, double T, double Tdot)
{
  double dg;
  flow_->y(s.data(), alpha.rawptr(), T, dg);

  Symmetric dkap = dkappa(edot, T);

  History hr = gather_blank_history_();
  flow_->h(s.data(), alpha.rawptr(), T, hr.rawptr());
  hr *= dg;
  
  History d_adot = gather_blank_history_().derivative<Symmetric>();
  d_adot.zero();
  outer_vec(hr.rawptr(), nhist(), dkap.data(), 6, d_adot.rawptr());

  flow_->h_time(s.data(), alpha.rawptr(), T, hr.rawptr());
  outer_update(hr.rawptr(), nhist(), dkap.data(), 6, d_adot.rawptr());

  return d_adot;
}

double WalkerKremplSwitchRule::work_rate(const Symmetric & s, const History & alpha, 
                                         const Symmetric & edot, double T, double Tdot)
{
  Symmetric erate;

  double kap = kappa(edot, T);

  Symmetric dir;
  double yv;
  flow_->g(s.data(), alpha.rawptr(), T, dir.s());
  flow_->y(s.data(), alpha.rawptr(), T, yv);

  erate += yv * kap * dir;

  flow_->g_temp(s.data(), alpha.rawptr(), T, dir.s());
  erate += Tdot * dir;

  flow_->g_time(s.data(), alpha.rawptr(), T, dir.s());
  erate += dir * kap;
  
  return s.contract(erate);
}

void WalkerKremplSwitchRule::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
}

double WalkerKremplSwitchRule::kappa(const Symmetric & edot, double T)
{
  double de = std::sqrt(2.0/3.0) * edot.dev().norm();
  return 1.0 - lambda_->value(T) + lambda_->value(T) * de / eps0_;
}

Symmetric WalkerKremplSwitchRule::dkappa(const Symmetric & edot, double T)
{
  Symmetric edev = edot.dev();
  double nde = edev.norm();

  if (nde == 0.0) {
    return Symmetric();
  }

  double fact = lambda_->value(T)  / eps0_ * std::sqrt(2.0/3.0) / nde;

  return edev * fact;
}

void WalkerKremplSwitchRule::override_guess(double * const x)
{
  flow_->override_guess(x);
}

SofteningModel::SofteningModel(ParameterSet & params) :
  NEMLObject(params)
{

}

std::string SofteningModel::type()
{
  return "SofteningModel";
}

std::unique_ptr<NEMLObject> SofteningModel::initialize(ParameterSet & params)
{
  return neml::make_unique<SofteningModel>(params); 
}

ParameterSet SofteningModel::parameters()
{
  ParameterSet pset(SofteningModel::type());

  return pset;
}

double SofteningModel::phi(double alpha, double T) const
{
  return 1.0;
}

double SofteningModel::dphi(double alpha, double T) const
{
  return 0.0;
}

WalkerSofteningModel::WalkerSofteningModel(ParameterSet & params) :
    SofteningModel(params),
    phi_0_(params.get_object_parameter<Interpolate>("phi_0")),
    phi_1_(params.get_object_parameter<Interpolate>("phi_1"))
{

}

std::string WalkerSofteningModel::type()
{
  return "WalkerSofteningModel";
}

std::unique_ptr<NEMLObject> WalkerSofteningModel::initialize(ParameterSet & params)
{
  return neml::make_unique<WalkerSofteningModel>(params); 
}

ParameterSet WalkerSofteningModel::parameters()
{
  ParameterSet pset(WalkerSofteningModel::type());

  pset.add_parameter<NEMLObject>("phi_0");
  pset.add_parameter<NEMLObject>("phi_1");

  return pset;
}

double WalkerSofteningModel::phi(double alpha, double T) const
{
  if (alpha <= 0.0) 
    return 1.0;
  else if (alpha < ainc_) 
    return phi_0_->value(T) * std::pow(ainc_, phi_1_->value(T)) / ainc_ * alpha +
        1.0;
  else 
    return 1.0 + phi_0_->value(T) * std::pow(alpha, phi_1_->value(T));
  
}

double WalkerSofteningModel::dphi(double alpha, double T) const
{
  if (alpha <= 0.0) 
   return phi_0_->value(T) * std::pow(ainc_, phi_1_->value(T)) / ainc_;
  else if (alpha < ainc_) 
    return phi_0_->value(T) * std::pow(ainc_, phi_1_->value(T)) / ainc_;
  else 
    return phi_1_->value(T) * phi_0_->value(T) * std::pow(alpha, phi_1_->value(T) - 1.0);
}

ThermalScaling::ThermalScaling(ParameterSet & params) :
    NEMLObject(params)
{

}

std::string ThermalScaling::type()
{
  return "ThermalScaling";
}

std::unique_ptr<NEMLObject> ThermalScaling::initialize(ParameterSet & params)
{
  return neml::make_unique<ThermalScaling>(params); 
}

ParameterSet ThermalScaling::parameters()
{
  ParameterSet pset(ThermalScaling::type());

  return pset;
}

double ThermalScaling::value(double T) const
{
  return 1.0;
}

ArrheniusThermalScaling::ArrheniusThermalScaling(ParameterSet & params) :
    ThermalScaling(params),
    Q_(params.get_object_parameter<Interpolate>("Q")),
    R_(params.get_parameter<double>("R")),
    T_ref_(params.get_parameter<double>("T_ref"))
{

}

std::string ArrheniusThermalScaling::type()
{
  return "ArrheniusThermalScaling";
}

std::unique_ptr<NEMLObject> ArrheniusThermalScaling::initialize(
    ParameterSet & params)
{
  return neml::make_unique<ArrheniusThermalScaling>(params); 
}

ParameterSet ArrheniusThermalScaling::parameters()
{
  ParameterSet pset(ArrheniusThermalScaling::type());

  pset.add_parameter<NEMLObject>("Q");
  pset.add_parameter<double>("R");
  pset.add_parameter<double>("T_ref");

  return pset;
}

double ArrheniusThermalScaling::value(double T) const
{
  return arr_(T) / arr_(T_ref_);
}

double ArrheniusThermalScaling::arr_(double T) const
{
  return std::exp(-Q_->value(T) / (R_ * T));
}

IsotropicHardening::IsotropicHardening(ParameterSet & params) :
    ScalarInternalVariable(params), 
    scale_(params.get_object_parameter<ThermalScaling>("scaling"))
{

}

/// Return zero for time rate by default 
double IsotropicHardening::ratet(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double IsotropicHardening::d_ratet_d_h(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double IsotropicHardening::d_ratet_d_a(VariableState & state) 
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double IsotropicHardening::d_ratet_d_adot(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double IsotropicHardening::d_ratet_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
Symmetric IsotropicHardening::d_ratet_d_s(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric IsotropicHardening::d_ratet_d_g(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for temperature rate by default 
double IsotropicHardening::rateT(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double IsotropicHardening::d_rateT_d_h(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double IsotropicHardening::d_rateT_d_a(VariableState & state) 
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double IsotropicHardening::d_rateT_d_adot(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double IsotropicHardening::d_rateT_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
Symmetric IsotropicHardening::d_rateT_d_s(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric IsotropicHardening::d_rateT_d_g(VariableState & state)
{
  return Symmetric::zero();
}


ConstantIsotropicHardening::ConstantIsotropicHardening(ParameterSet & params) :
      IsotropicHardening(params)
{

}

std::string ConstantIsotropicHardening::type()
{
  return "ConstantIsotropicHardening";
}

ParameterSet ConstantIsotropicHardening::parameters()
{
  ParameterSet pset(ConstantIsotropicHardening::type());
  
  pset.add_optional_parameter<std::string>("name", "R");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          default_scaling());

  return pset;
}

std::unique_ptr<NEMLObject> ConstantIsotropicHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<ConstantIsotropicHardening>(params); 
}

double ConstantIsotropicHardening::initial_value()
{
  return 0;
}

double ConstantIsotropicHardening::ratep(VariableState & state)
{
  return 0;
}

double ConstantIsotropicHardening::d_ratep_d_h(VariableState & state)
{
  return 0;
}

double ConstantIsotropicHardening::d_ratep_d_a(VariableState & state)
{
  return 0;
}

double ConstantIsotropicHardening::d_ratep_d_adot(VariableState & state)
{
  return 0;
}

double ConstantIsotropicHardening::d_ratep_d_D(VariableState & state)
{
  return 0;
}

Symmetric ConstantIsotropicHardening::d_ratep_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric ConstantIsotropicHardening::d_ratep_d_g(VariableState & state)
{
  return Symmetric();
}


WalkerIsotropicHardening::WalkerIsotropicHardening(ParameterSet & params) :
      IsotropicHardening(params),
      r0_(params.get_object_parameter<Interpolate>("r0")),
      Rinf_(params.get_object_parameter<Interpolate>("Rinf")),
      R0_(params.get_object_parameter<Interpolate>("R0")),
      r1_(params.get_object_parameter<Interpolate>("r1")), 
      r2_(params.get_object_parameter<Interpolate>("r2"))
{

}

std::string WalkerIsotropicHardening::type()
{
  return "WalkerIsotropicHardening";
}

ParameterSet WalkerIsotropicHardening::parameters()
{
  ParameterSet pset(WalkerIsotropicHardening::type());

  pset.add_parameter<NEMLObject>("r0");
  pset.add_parameter<NEMLObject>("Rinf");
  pset.add_parameter<NEMLObject>("R0");
  pset.add_parameter<NEMLObject>("r1");
  pset.add_parameter<NEMLObject>("r2");
  pset.add_optional_parameter<std::string>("name", "R");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          default_scaling());

  return pset;
}

std::unique_ptr<NEMLObject> WalkerIsotropicHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<WalkerIsotropicHardening>(params); 
}

double WalkerIsotropicHardening::initial_value()
{
  return 0;
}

double WalkerIsotropicHardening::ratep(VariableState & state)
{
  return r0_->value(state.T) * (Rinf_->value(state.T) - state.h);
}

double WalkerIsotropicHardening::d_ratep_d_h(VariableState & state)
{
  return -r0_->value(state.T);
}

double WalkerIsotropicHardening::d_ratep_d_a(VariableState & state)
{
  return 0;
}

double WalkerIsotropicHardening::d_ratep_d_adot(VariableState & state)
{
  return 0;
}

double WalkerIsotropicHardening::d_ratep_d_D(VariableState & state)
{
  return 0;
}

Symmetric WalkerIsotropicHardening::d_ratep_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric WalkerIsotropicHardening::d_ratep_d_g(VariableState & state)
{
  return Symmetric();
}

double WalkerIsotropicHardening::ratet(VariableState & state)
{
  double d = R0_->value(state.T) - state.h;
  return r1_->value(state.T) * d * std::pow(std::fabs(d), r2_->value(state.T) - 1.0);
}

double WalkerIsotropicHardening::d_ratet_d_h(VariableState & state)
{
  return -r1_->value(state.T) * r2_->value(state.T) * 
      std::pow(std::fabs(R0_->value(state.T) - state.h), r2_->value(state.T) - 1.0);
}

double WalkerIsotropicHardening::d_ratet_d_a(VariableState & state)
{
  return 0;
}

double WalkerIsotropicHardening::d_ratet_d_adot(VariableState & state)
{
  return 0;
}

double WalkerIsotropicHardening::d_ratet_d_D(VariableState & state)
{
  return 0;
}

Symmetric WalkerIsotropicHardening::d_ratet_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric WalkerIsotropicHardening::d_ratet_d_g(VariableState & state)
{
  return Symmetric();
}

DragStress::DragStress(ParameterSet & params) :
    ScalarInternalVariable(params),
    scale_(params.get_object_parameter<ThermalScaling>("scaling"))
{

}

/// Makes no sense in this context
double DragStress::d_ratep_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for time rate by default 
double DragStress::ratet(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double DragStress::d_ratet_d_h(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double DragStress::d_ratet_d_a(VariableState & state) 
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double DragStress::d_ratet_d_adot(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
double DragStress::d_ratet_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for the time rate derivatives by default
Symmetric DragStress::d_ratet_d_s(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric DragStress::d_ratet_d_g(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for temperature rate by default 
double DragStress::rateT(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double DragStress::d_rateT_d_h(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double DragStress::d_rateT_d_a(VariableState & state) 
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double DragStress::d_rateT_d_adot(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
double DragStress::d_rateT_d_D(VariableState & state)
{
  return 0;
}

/// Return zero for the temperature rate derivatives by default
Symmetric DragStress::d_rateT_d_s(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric DragStress::d_rateT_d_g(VariableState & state)
{
  return Symmetric::zero();
}

ConstantDragStress::ConstantDragStress(ParameterSet & params) :
      DragStress(params), 
      value_(params.get_parameter<double>("value"))
{

}

std::string ConstantDragStress::type()
{
  return "ConstantDragStress";
}

ParameterSet ConstantDragStress::parameters()
{
  ParameterSet pset(ConstantDragStress::type());
  
  pset.add_parameter<double>("value");
  pset.add_optional_parameter<std::string>("name", "D");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          default_scaling());

  return pset;
}

std::unique_ptr<NEMLObject> ConstantDragStress::initialize(
    ParameterSet & params)
{
  return neml::make_unique<ConstantDragStress>(params); 
}

double ConstantDragStress::initial_value()
{
  return value_;
}

double ConstantDragStress::D_xi(double T)
{
  return 1; // Can be an arbitrary value for this model
}

double ConstantDragStress::D_0(double T)
{
  return value_; // Straightforward!
}

double ConstantDragStress::ratep(VariableState & state)
{
  return 0;
}

double ConstantDragStress::d_ratep_d_h(VariableState & state)
{
  return 0;
}

double ConstantDragStress::d_ratep_d_a(VariableState & state)
{
  return 0;
}

double ConstantDragStress::d_ratep_d_adot(VariableState & state)
{
  return 0;
}

Symmetric ConstantDragStress::d_ratep_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric ConstantDragStress::d_ratep_d_g(VariableState & state)
{
  return Symmetric();
}

WalkerDragStress::WalkerDragStress(ParameterSet & params) :
      DragStress(params),
      d0_(params.get_object_parameter<Interpolate>("d0")), 
      d1_(params.get_object_parameter<Interpolate>("d1")), 
      d2_(params.get_object_parameter<Interpolate>("d2")),
      D_xi_(params.get_object_parameter<Interpolate>("D_xi")), 
      D_0_(params.get_parameter<double>("D_0")), 
      softening_(params.get_object_parameter<SofteningModel>("softening"))
{

}

std::string WalkerDragStress::type()
{
  return "WalkerDragStress";
}

ParameterSet WalkerDragStress::parameters()
{
  ParameterSet pset(WalkerDragStress::type());
  
  pset.add_parameter<NEMLObject>("d0");
  pset.add_parameter<NEMLObject>("d1");
  pset.add_parameter<NEMLObject>("d2");
  pset.add_parameter<NEMLObject>("D_xi");
  pset.add_parameter<double>("D_0");
  pset.add_parameter<NEMLObject>("softening");
  pset.add_optional_parameter<std::string>("name", "D");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          default_scaling());

  return pset;
}

std::unique_ptr<NEMLObject> WalkerDragStress::initialize(
    ParameterSet & params)
{
  return neml::make_unique<WalkerDragStress>(params); 
}

double WalkerDragStress::initial_value()
{
  return D_0_;
}

double WalkerDragStress::D_xi(double T)
{
  return D_xi_->value(T); // Explicit parameter
}

double WalkerDragStress::D_0(double T)
{
  return D_0_; // Straightforward!
}

double WalkerDragStress::ratep(VariableState & state)
{
  return d0_->value(state.T) * (1.0 - (state.h - D_0_) / D_xi_->value(state.T));
}

double WalkerDragStress::d_ratep_d_h(VariableState & state)
{
  return -d0_->value(state.T) / D_xi_->value(state.T);
}

double WalkerDragStress::d_ratep_d_a(VariableState & state)
{
  return 0;
}

double WalkerDragStress::d_ratep_d_adot(VariableState & state)
{
  return 0;
}

Symmetric WalkerDragStress::d_ratep_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric WalkerDragStress::d_ratep_d_g(VariableState & state)
{
  return Symmetric();
}

double WalkerDragStress::ratet(VariableState & state)
{
  if ((state.h-D_0_) <= 0.0)
    return 0;
  else
    return -scale_->value(state.T) * softening_->phi(state.a, state.T) * 
        d1_->value(state.T) * std::pow(state.h - D_0_, d2_->value(state.T));
}

double WalkerDragStress::d_ratet_d_h(VariableState & state)
{
  if ((state.h-D_0_) <= 0.0)
    return 1;
  else
    return -d2_->value(state.T) * scale_->value(state.T) * softening_->phi(state.a, state.T) * 
        d1_->value(state.T) * std::pow(state.h - D_0_, d2_->value(state.T) - 1.0);
}

double WalkerDragStress::d_ratet_d_a(VariableState & state)
{
  if ((state.h-D_0_) <= 0.0)
    return 1;
  else
    return -scale_->value(state.T) * softening_->dphi(state.a, state.T) * 
        d1_->value(state.T) * std::pow(state.h - D_0_, d2_->value(state.T));
}

double WalkerDragStress::d_ratet_d_adot(VariableState & state)
{
  return 0;
}

Symmetric WalkerDragStress::d_ratet_d_s(VariableState & state)
{
  return Symmetric();
}

Symmetric WalkerDragStress::d_ratet_d_g(VariableState & state)
{
  return Symmetric();
}

KinematicHardening::KinematicHardening(ParameterSet & params) :
    SymmetricInternalVariable(params), 
    scale_(params.get_object_parameter<ThermalScaling>("scaling"))
{

}

/// Return zero for time rate by default 
Symmetric KinematicHardening::ratet(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
SymSymR4 KinematicHardening::d_ratet_d_h(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric KinematicHardening::d_ratet_d_a(VariableState & state) 
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric KinematicHardening::d_ratet_d_adot(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
Symmetric KinematicHardening::d_ratet_d_D(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the time rate derivatives by default
SymSymR4 KinematicHardening::d_ratet_d_s(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for the time rate derivatives by default
SymSymR4 KinematicHardening::d_ratet_d_g(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for temperature rate by default 
Symmetric KinematicHardening::rateT(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
SymSymR4 KinematicHardening::d_rateT_d_h(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric KinematicHardening::d_rateT_d_a(VariableState & state) 
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric KinematicHardening::d_rateT_d_adot(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
Symmetric KinematicHardening::d_rateT_d_D(VariableState & state)
{
  return Symmetric::zero();
}

/// Return zero for the temperature rate derivatives by default
SymSymR4 KinematicHardening::d_rateT_d_s(VariableState & state)
{
  return SymSymR4::zero();
}

/// Return zero for the temperature rate derivatives by default
SymSymR4 KinematicHardening::d_rateT_d_g(VariableState & state)
{
  return SymSymR4::zero();
}

FAKinematicHardening::FAKinematicHardening(ParameterSet & params) : 
    KinematicHardening(params),
    c_(params.get_object_parameter<Interpolate>("c")),
    g_(params.get_object_parameter<Interpolate>("g"))
{

}

std::string FAKinematicHardening::type()
{
  return "FAKinematicHardening";
}

ParameterSet FAKinematicHardening::parameters()
{
  ParameterSet pset(FAKinematicHardening::type());
  
  pset.add_parameter<NEMLObject>("c");
  pset.add_parameter<NEMLObject>("g");
  pset.add_optional_parameter<std::string>("name", "X");
  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          default_scaling());

  return pset;
}

std::unique_ptr<NEMLObject> FAKinematicHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<FAKinematicHardening>(params); 
}

Symmetric FAKinematicHardening::initial_value()
{
  return Symmetric::zero();
}

Symmetric FAKinematicHardening::ratep(VariableState & state)
{
  return 2.0/3.0 * c_->value(state.T) * state.g - g_->value(state.T) * state.h;
}

SymSymR4 FAKinematicHardening::d_ratep_d_h(VariableState & state)
{
  return -g_->value(state.T) * SymSymR4::id();
}

Symmetric FAKinematicHardening::d_ratep_d_a(VariableState & state) 
{
  return Symmetric::zero();
}

Symmetric FAKinematicHardening::d_ratep_d_adot(VariableState & state)
{
  return Symmetric::zero();
}

Symmetric FAKinematicHardening::d_ratep_d_D(VariableState & state)
{
  return Symmetric::zero();
}

SymSymR4 FAKinematicHardening::d_ratep_d_s(VariableState & state)
{
  return SymSymR4::zero();
}

SymSymR4 FAKinematicHardening::d_ratep_d_g(VariableState & state)
{
  return 2.0/3.0 * c_->value(state.T) * SymSymR4::id();
}


WalkerKinematicHardening::WalkerKinematicHardening(ParameterSet & params) : 
      KinematicHardening(params), 
      c0_(params.get_object_parameter<Interpolate>("c0")),
      c1_(params.get_object_parameter<Interpolate>("c1")), 
      c2_(params.get_object_parameter<Interpolate>("c2")), 
      l0_(params.get_object_parameter<Interpolate>("l0")), 
      l1_(params.get_object_parameter<Interpolate>("l1")),
      l_(params.get_object_parameter<Interpolate>("l")), 
      b0_(params.get_object_parameter<Interpolate>("b0")), 
      x0_(params.get_object_parameter<Interpolate>("x0")), 
      x1_(params.get_object_parameter<Interpolate>("x1")), 
      softening_(params.get_object_parameter<SofteningModel>("softening"))
{

}

std::string WalkerKinematicHardening::type()
{
  return "WalkerKinematicHardening";
}

ParameterSet WalkerKinematicHardening::parameters()
{
  ParameterSet pset(WalkerKinematicHardening::type());
  
  pset.add_parameter<NEMLObject>("c0");
  pset.add_parameter<NEMLObject>("c1");
  pset.add_parameter<NEMLObject>("c2");
  pset.add_parameter<NEMLObject>("l0");
  pset.add_parameter<NEMLObject>("l1");
  pset.add_parameter<NEMLObject>("l");
  pset.add_parameter<NEMLObject>("b0");
  pset.add_parameter<NEMLObject>("x0");
  pset.add_parameter<NEMLObject>("x1");
  pset.add_parameter<NEMLObject>("softening");

  pset.add_optional_parameter<NEMLObject>("scaling", 
                                          default_scaling());
  pset.add_optional_parameter<std::string>("name", "X");

  return pset;
}

std::unique_ptr<NEMLObject> WalkerKinematicHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<WalkerKinematicHardening>(params); 
}

Symmetric WalkerKinematicHardening::initial_value()
{
  return Symmetric::zero();
}

Symmetric WalkerKinematicHardening::ratep(VariableState & state)
{
  return c_(state) * (2.0/3.0 * state.g - b_(state) / L_(state));
}

SymSymR4 WalkerKinematicHardening::d_ratep_d_h(VariableState & state)
{
  return -c_(state) * db_dx_(state) / L_(state);
}

Symmetric WalkerKinematicHardening::d_ratep_d_a(VariableState & state) 
{
  return c_(state) * b_(state) / std::pow(L_(state), 2.0) * dL_(state);
}

Symmetric WalkerKinematicHardening::d_ratep_d_adot(VariableState & state)
{
  return dc_(state) * (2.0/3.0 * state.g - b_(state) / L_(state));;
}

Symmetric WalkerKinematicHardening::d_ratep_d_D(VariableState & state)
{
  return Symmetric::zero();
}

SymSymR4 WalkerKinematicHardening::d_ratep_d_s(VariableState & state)
{
  return -c_(state) * db_ds_(state) / L_(state);
}

SymSymR4 WalkerKinematicHardening::d_ratep_d_g(VariableState & state)
{
  return c_(state) * 2.0/3.0 * SymSymR4::id();
}

Symmetric WalkerKinematicHardening::ratet(VariableState & state)
{
  if ((state.h.norm() == 0.0) || (state.D <= 0.0)) 
    return Symmetric::zero();
  else 
    return -scale_->value(state.T) * x0_->value(state.T) * 
        softening_->phi(state.a, state.T) * 
        std::pow(std::sqrt(3.0/2.0) * state.h.norm() / state.D, x1_->value(state.T))
        * state.h / (state.h.norm() * std::sqrt(3.0/2));
}

SymSymR4 WalkerKinematicHardening::d_ratet_d_h(VariableState & state)
{
  if ((state.h.norm() == 0) || (state.D <= 0.0)) {
    return SymSymR4::zero();
  }
  else {
    Symmetric d = state.h / state.h.norm();
    SymSymR4 dd = 1.0 / (std::sqrt(3.0/2.0) * state.h.norm()) * (SymSymR4::id() - 
                                                                 douter(d, d));
    return -scale_->value(state.T) * x0_->value(state.T) * 
        softening_->phi(state.a, state.T) * 
        (x1_->value(state.T)/state.D * std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D,
                                 x1_->value(state.T) - 1.0)  * douter(d,d)
         + std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D,
                    x1_->value(state.T)) * dd);
  }
}

Symmetric WalkerKinematicHardening::d_ratet_d_a(VariableState & state) 
{
  if ((state.h.norm() == 0) || (state.D <= 0.0))
    return Symmetric::zero();
  else
    return -scale_->value(state.T) * x0_->value(state.T) * 
        softening_->dphi(state.a, state.T) * 
       std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D, x1_->value(state.T))
       * state.h / (state.h.norm() * std::sqrt(3.0/2));
}

Symmetric WalkerKinematicHardening::d_ratet_d_adot(VariableState & state)
{
  return Symmetric::zero();
}

Symmetric WalkerKinematicHardening::d_ratet_d_D(VariableState & state)
{
  if ((state.h.norm() == 0) || (state.D <= 0.0))
    return Symmetric::zero();
  else
    return scale_->value(state.T) * x0_->value(state.T) * 
        softening_->phi(state.a, state.T) * x1_->value(state.T) *
        std::pow(std::sqrt(3.0/2) * state.h.norm() / state.D,
                 x1_->value(state.T)-1.0)
        * state.h / (state.h.norm() * std::sqrt(3.0/2))
        * std::sqrt(3.0/2) * state.h.norm() / (state.D*state.D);
}

SymSymR4 WalkerKinematicHardening::d_ratet_d_s(VariableState & state)
{
  return SymSymR4::zero();
}

SymSymR4 WalkerKinematicHardening::d_ratet_d_g(VariableState & state)
{
  return SymSymR4::zero();
}

double WalkerKinematicHardening::c_(VariableState & state)
{
  if (state.adot <= 0.0)
    return c0_->value(state.T);
  else
    return c0_->value(state.T) + c1_->value(state.T) * 
        std::pow(state.adot, 1.0/c2_->value(state.T));
}

double WalkerKinematicHardening::dc_(VariableState & state)
{
  if (state.adot <= 0.0)
    return 0.0;
  else
    return c1_->value(state.T)/c2_->value(state.T) * 
        std::pow(state.adot, 1.0/c2_->value(state.T)-1.0);
}

double WalkerKinematicHardening::L_(VariableState & state)
{
  if (state.a <= 0.0) return l1_->value(state.T) + 1.0;
  else return l_->value(state.T) * (l1_->value(state.T) + (1.0 - l1_->value(state.T))
                               * std::exp(-l0_->value(state.T) * state.a));
}

double WalkerKinematicHardening::dL_(VariableState & state)
{
  if (state.a <= 0.0) return -l0_->value(state.T) * l_->value(state.T);
  else return -l0_->value(state.T) * l_->value(state.T) * (1.0 - l1_->value(state.T))
                                * std::exp(-l0_->value(state.T) * state.a);
}

Symmetric WalkerKinematicHardening::n_(VariableState & state)
{
  if ((state.s.dev() - state.h).norm() == 0.0) 
    return Symmetric::zero();
  else
    return (3.0/2.0) * (state.s.dev() - state.h) / 
        (std::sqrt(3.0/2.0) * (state.s.dev() - state.h).norm());
}

Symmetric WalkerKinematicHardening::b_(VariableState & state)
{
  Symmetric n = n_(state);
  return (1.0 - b0_->value(state.T)) * state.h + 2.0/3.0 * b0_->value(state.T) *
      douter(n,n).dot(state.h);
}

SymSymR4 WalkerKinematicHardening::db_ds_(VariableState & state)
{
  SymSymR4 Nbar = dN_(state).dot(SymSymR4::id_dev());
  Symmetric n = n_(state);

  return 2.0/3.0 * b0_->value(state.T) * (Nbar * n.contract(state.h) + 
                                          douter(n,
                                                 Nbar.dot(state.h).transpose()));
}

SymSymR4 WalkerKinematicHardening::db_dx_(VariableState & state)
{
  Symmetric n = n_(state);
  SymSymR4 N = dN_(state);
  return (1.0 - b0_->value(state.T)) * SymSymR4::id() 
      + 2.0/3.0 * b0_->value(state.T) * douter(n,n)
      - 2.0/3.0 * b0_->value(state.T) * (N * n.contract(state.h) + 
                                         douter(n, 
                                                N.dot(state.h).transpose()));
}

SymSymR4 WalkerKinematicHardening::dN_(VariableState & state)
{
  Symmetric d = state.s.dev() - state.h;
  double dn = d.norm();
  if (dn == 0.0)
    return SymSymR4::id();
  else
    return std::sqrt(3.0/2.0) / dn * (SymSymR4::id() - douter(d/dn,d/dn));
}

WrappedViscoPlasticFlowRule::WrappedViscoPlasticFlowRule(ParameterSet & params) :
    ViscoPlasticFlowRule(params),
    stored_hist_(false)
{

}

History WrappedViscoPlasticFlowRule::blank_hist_() const
{
  return stored_hist_;
}

History WrappedViscoPlasticFlowRule::create_blank_hist_() const
{
  History h;
  populate_hist(h);
  return h;
}

History WrappedViscoPlasticFlowRule::gather_hist_(double * const h) const
{
  History hv = blank_hist_();
  hv.set_data(h);
  return hv;
}

History WrappedViscoPlasticFlowRule::gather_hist_(const double * const h) const
{
  History hv = blank_hist_();
  hv.set_data(const_cast<double*>(h));
  return hv;
}

State WrappedViscoPlasticFlowRule::make_state_(const double * const s, const double *
                                               const alpha, double T) const
{
  return State(Symmetric(s), gather_hist_(alpha), T);
}

// Rate rule
void WrappedViscoPlasticFlowRule::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  y(make_state_(s, alpha, T), yv);
}

void WrappedViscoPlasticFlowRule::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  Symmetric res(dyv);
  dy_ds(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  // This is transposed, but that doesn't matter as it's flat
  History res = gather_derivative_<double>(dyv);
  dy_da(make_state_(s, alpha, T), res);
}

// Flow rule
void WrappedViscoPlasticFlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  Symmetric res(gv);
  g(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  SymSymR4 res(dgv);
  dg_ds(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::dg_da(const double * const s, const double * const alpha, double T,
             double * const dgv) const
{
  // This is transposed and it does matter
  double * temp = new double [nhist() * 6];
  History res = gather_derivative_<Symmetric>(temp);
  dg_da(make_state_(s, alpha, T), res);
  
  for (size_t i = 0; i < nhist(); i++)
    for (size_t j = 0; j < 6; j++)
      dgv[CINDEX(j,i,nhist())] = temp[CINDEX(i,j,6)];
  delete [] temp;
}

// Hardening rule
void WrappedViscoPlasticFlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  History res = gather_hist_(hv);
  h(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<Symmetric>(dhv);
  dh_ds(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  // Very messed up
  double * temp = new double[nhist() * nhist()];
  History res = gather_derivative_<History>(temp);

  dh_da(make_state_(s, alpha, T), res);

  res.unravel_hh(blank_hist_(), dhv);

  delete [] temp;
}

// Hardening rule wrt time
void WrappedViscoPlasticFlowRule::h_time(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  History res = gather_hist_(hv);
  h_time(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::h_time(const State & state, History & res) const
{
  res.zero();
}

void WrappedViscoPlasticFlowRule::dh_ds_time(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<Symmetric>(dhv);
  dh_ds_time(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::dh_ds_time(const State & state, History & res) const
{
  res.zero();
}

void WrappedViscoPlasticFlowRule::dh_da_time(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  // Very messed up
  double * temp = new double[nhist() * nhist()];
  History res = gather_derivative_<History>(temp);

  dh_da_time(make_state_(s, alpha, T), res);

  res.unravel_hh(blank_hist_(), dhv);

  delete [] temp;
}

void WrappedViscoPlasticFlowRule::dh_da_time(const State & state, History & res) const
{
  res.zero();
}

// Hardening rule wrt temperature
void WrappedViscoPlasticFlowRule::h_temp(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  History res = gather_hist_(hv);
  h_temp(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::h_temp(const State & state, History & res) const
{
  res.zero();
}

void WrappedViscoPlasticFlowRule::dh_ds_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<Symmetric>(dhv);
  dh_ds_temp(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::dh_ds_temp(const State & state, History & res) const
{
  res.zero();
}

void WrappedViscoPlasticFlowRule::dh_da_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  History res = gather_derivative_<History>(dhv);
  dh_da_temp(make_state_(s, alpha, T), res);
}

void WrappedViscoPlasticFlowRule::dh_da_temp(const State & state, History & res) const
{
  res.zero();
}

TestFlowRule::TestFlowRule(ParameterSet & params)
  : WrappedViscoPlasticFlowRule(params), 
    eps0_(params.get_parameter<double>("eps0")), 
    D_(params.get_parameter<double>(prefix("D"))), 
    n_(params.get_parameter<double>("n")), 
    s0_(params.get_parameter<double>("s0")), 
    K_(params.get_parameter<double>("K"))
{
  populate_hist(stored_hist_);
}

std::string TestFlowRule::type()
{
  return "TestFlowRule";
}

std::unique_ptr<NEMLObject> TestFlowRule::initialize(ParameterSet & params)
{
  return neml::make_unique<TestFlowRule>(params);
}

ParameterSet TestFlowRule::parameters()
{
  ParameterSet pset(TestFlowRule::type());

  pset.add_parameter<double>("eps0");
  pset.add_parameter<double>("D");
  pset.add_parameter<double>("n");
  pset.add_parameter<double>("s0");
  pset.add_parameter<double>("K");

  return pset;
}

void TestFlowRule::populate_hist(History & h) const
{
  h.add<double>(prefix("alpha"));
  h.add<double>(prefix("iso"));
}

void TestFlowRule::init_hist(History & h) const
{
  h.get<double>(prefix("alpha")) = 0.0;
  h.get<double>(prefix("iso")) = s0_;
}

void TestFlowRule::y(const State & state, double & res) const
{
  double h = (std::sqrt(3.0/2.0) * state.S.dev().norm() - state.h.get<double>(prefix("iso")))
      / D_;
  if (h > 0.0) {
    res = eps0_ * std::pow(h, n_);
  }
  else {
    res = 0;
  }
}

void TestFlowRule::dy_ds(const State & state, Symmetric & res) const
{
  double h = (std::sqrt(3.0/2.0) * state.S.dev().norm() - state.h.get<double>(prefix("iso")))
      / D_;
  if (h > 0.0) {
    res = eps0_ * n_ * std::pow(h, n_-1.0) * std::sqrt(3.0/2.0) * state.S.dev()
        / state.S.dev().norm() / D_;
  }
  else {
    res = Symmetric::zero();
  }
}

void TestFlowRule::dy_da(const State & state, History & res) const
{
  double h = (std::sqrt(3.0/2.0) * state.S.dev().norm() - state.h.get<double>(prefix("iso")))
      / D_;
  res.zero();
  if (h > 0.0) {
    res.get<double>(prefix("iso")) = -eps0_ * n_ * std::pow(h, n_-1.0) / D_;
  }
}

void TestFlowRule::g(const State & state, Symmetric & res) const
{
  double ns = state.S.dev().norm();
  if (ns > 0)
    res = std::sqrt(3.0 / 2.0) * state.S.dev() / ns;
  else
    res = Symmetric::zero();
}

void TestFlowRule::dg_ds(const State & state, SymSymR4 & res) const
{
  Symmetric s = state.S.dev();
  double ns = s.norm();
  Symmetric sn = s / ns;
  res = std::sqrt(3.0 / 2.0) / ns * (SymSymR4::id_dev() - douter(sn, sn));
}

void TestFlowRule::dg_da(const State & state, History & res) const
{
  res.zero();
}

void TestFlowRule::h(const State & state, History & res) const
{
  res.get<double>(prefix("alpha")) = 1.0;
  res.get<double>(prefix("iso")) = K_;
}

void TestFlowRule::dh_ds(const State & state, History & res) const
{
  res.zero();
}

void TestFlowRule::dh_da(const State & state, History & res) const
{
  res.zero();
}

WalkerFlowRule::WalkerFlowRule(ParameterSet & params) :
    WrappedViscoPlasticFlowRule(params),
    eps0_(params.get_object_parameter<Interpolate>("eps0")),
    softening_(params.get_object_parameter<SofteningModel>("softening")),
    scaling_(params.get_object_parameter<ThermalScaling>("scaling")),
    n_(params.get_object_parameter<Interpolate>("n")),
    k_(params.get_object_parameter<Interpolate>("k")),
    m_(params.get_object_parameter<Interpolate>("m")),
    R_(params.get_object_parameter<IsotropicHardening>(prefix("R"))),
    D_(params.get_object_parameter<DragStress>(prefix("D"))),
    X_(params.get_object_parameter_vector<KinematicHardening>("X"))
{
  // The variable names to something canonical
  R_->set_name(prefix("R"));
  D_->set_name(prefix("D"));
  int i = 0;
  for (auto X : X_) {
    X->set_name("X" + std::to_string(i));
    i++;
  }

  // Set the thermal scaling models to be the same
  R_->set_scaling(scaling_);
  D_->set_scaling(scaling_);
  for (auto X : X_) {
    X->set_scaling(scaling_);
  }

  populate_hist(stored_hist_);
}

std::string WalkerFlowRule::type()
{
  return "WalkerFlowRule";
}

std::unique_ptr<NEMLObject> WalkerFlowRule::initialize(ParameterSet & params)
{
  return neml::make_unique<WalkerFlowRule>(params);
}

ParameterSet WalkerFlowRule::parameters()
{
  ParameterSet pset(WalkerFlowRule::type());

  pset.add_parameter<NEMLObject>("eps0");
  pset.add_parameter<NEMLObject>("softening");
  pset.add_parameter<NEMLObject>("scaling");
  pset.add_parameter<NEMLObject>("n");
  pset.add_parameter<NEMLObject>("k");
  pset.add_parameter<NEMLObject>("m");
  pset.add_parameter<NEMLObject>("R");
  pset.add_parameter<NEMLObject>("D");
  pset.add_parameter<std::vector<NEMLObject>>("X");

  return pset;
}

void WalkerFlowRule::populate_hist(History & h) const
{
  h.add<double>(prefix("alpha"));
  h.add<double>(prefix(R_->name()));
  h.add<double>(prefix(D_->name()));
  for (auto X : X_) {
    h.add<Symmetric>(prefix(X->name()));
  }
}

void WalkerFlowRule::init_hist(History & h) const
{
  h.get<double>(prefix("alpha")) = 0;
  h.get<double>(prefix(R_->name())) = R_->initial_value();
  h.get<double>(prefix(D_->name())) = D_->initial_value();
  for (auto X : X_) {
    h.get<Symmetric>(prefix(X->name())) = X->initial_value();
  }
}

void WalkerFlowRule::y(const State & state, double & res) const
{
  res = prefactor_(state) * flow_(state);
}

void WalkerFlowRule::dy_ds(const State & state, Symmetric & res) const
{
  Symmetric d = state.S.dev() - TX_(state);
  if (d.norm() == 0)
    res = Symmetric::zero();
  else
    res = prefactor_(state) * dflow_(state) * 
        std::sqrt(3.0/2.0)/(state.h.get<double>(prefix("D")) * d.norm()) * 
        SymSymR4::id_dev().dot(d);
}

void WalkerFlowRule::dy_da(const State & state, History & res) const
{
  // These are very annoying because it covers all the history variables...
  
  // alpha
  res.get<double>(prefix("alpha")) = eps0_->value(state.T) * 
      softening_->dphi(state.h.get<double>(prefix("alpha")), state.T) * 
      scaling_->value(state.T) * flow_(state);

  // R
  double D = state.h.get<double>(prefix("D"));
  double D0 = D_->D_0(state.T);
  double Dx = D_->D_xi(state.T);
  double m = m_->value(state.T);

  double Dr = 1.0e-8;
  if ((D - D0) > 0) Dr = (D-D0)/Dx;
  
  res.get<double>(prefix("R")) = -prefactor_(state) * dflow_(state) * 
      std::pow(Dr, m) / D;

  // D
  Symmetric d = state.S.dev() - TX_(state);
  double k = k_->value(state.T);
  double R = state.h.get<double>(prefix("R"));
  double p1 = -(k + R) * m * std::pow(Dr, m-1.0) / (Dx*D);
  double p2 = -(std::sqrt(3.0/2.0) * d.norm() - (k + R) * std::pow(Dr, m)) / (D*D);

  res.get<double>(prefix("D")) = prefactor_(state) * dflow_(state) * 
      (p1 + p2);

  // The backstresses
  for (auto X : X_)
    if (d.norm() == 0) {
      res.get<Symmetric>(prefix(X->name())) = Symmetric::zero();
    }
    else {
      res.get<Symmetric>(prefix(X->name())) = -prefactor_(state) * dflow_(state) * 
          std::sqrt(3.0/2.0)/(state.h.get<double>(prefix("D")) * d.norm()) * d;
    }
}

void WalkerFlowRule::g(const State & state, Symmetric & res) const
{
  Symmetric d = state.S.dev() - TX_(state);
  if (d.norm() == 0)
    res = Symmetric::zero();
  else
    res = 3.0/2.0 * d / (std::sqrt(3.0/2.0) * d.norm());
}

void WalkerFlowRule::dg_ds(const State & state, SymSymR4 & res) const
{
  SymSymR4 G = G_(state);
  res = G.dot(SymSymR4::id_dev());
}

void WalkerFlowRule::dg_da(const State & state, History & res) const
{
  res.get<Symmetric>(prefix("alpha")) = Symmetric::zero();
  res.get<Symmetric>(prefix("R")) = Symmetric::zero();
  res.get<Symmetric>(prefix("D")) = Symmetric::zero();

  SymSymR4 G = G_(state);
  for (auto X : X_)
    res.get<SymSymR4>(prefix(X->name())) = -G;
}

void WalkerFlowRule::h(const State & state, History & res) const
{
  // Scalar variables
  res.get<double>(prefix("alpha")) = 1.0;
  auto ss = scalar_state_(state);
  ss.h = state.h.get<double>(prefix("R"));
  res.get<double>(prefix("R")) = R_->ratep(ss);
  ss.h = state.h.get<double>(prefix("D"));
  res.get<double>(prefix("D")) = D_->ratep(ss);
  
  // Backstresses
  auto Ss = symmetric_state_(state);
  for (auto X : X_) {
    Ss.h.copy_data(state.h.get<Symmetric>(prefix(X->name())).data());
    res.get<Symmetric>(prefix(X->name())) = X->ratep(Ss);
  }
}

void WalkerFlowRule::dh_ds(const State & state, History & res) const
{
  // Again highly annoying, but at least there's a pattern here
  // Alpha
  res.get<Symmetric>(prefix("alpha")) = Symmetric::zero();
  
  // Common junk
  Symmetric dy;
  dy_ds(state, dy);
  SymSymR4 dg;
  dg_ds(state, dg);

  auto ss = scalar_state_(state);  

  // R
  ss.h = state.h.get<double>(prefix("R"));
  res.get<Symmetric>(prefix("R")) = 
        R_->d_ratep_d_s(ss) 
      + R_->d_ratep_d_adot(ss) * dy
      + dg.dot(R_->d_ratep_d_g(ss)).transpose();

  // D
  ss.h = state.h.get<double>(prefix("D"));
  res.get<Symmetric>(prefix("D")) = 
        D_->d_ratep_d_s(ss) 
      + D_->d_ratep_d_adot(ss) * dy
      + dg.dot(D_->d_ratep_d_g(ss)).transpose();

  // Backstresses
  auto Ss = symmetric_state_(state);
  for (auto X : X_) {
    Ss.h.copy_data(state.h.get<Symmetric>(prefix(X->name())).data());
    res.get<SymSymR4>(prefix(X->name())) = 
          X->d_ratep_d_s(Ss)
        + douter(X->d_ratep_d_adot(Ss), dy)
        + X->d_ratep_d_g(Ss).dot(dg);
  }
}

void WalkerFlowRule::dh_da(const State & state, History & res) const
{
  // Well this is the worst...
  
  // Common stuff we need...
  History dy = blank_derivative_<double>();
  dy_da(state, dy);

  History dg = blank_derivative_<Symmetric>();
  dg_da(state, dg);

  // Make sure the result starts out as zero
  res.zero();

  // Scalar variables
  auto ss = scalar_state_(state);

  // Alpha
  // (zero)

  // R
  ss.h = state.h.get<double>(prefix("R"));
  // a
  res.get<double>(dprefix("R", "alpha")) =
        R_->d_ratep_d_a(ss)
      + R_->d_ratep_d_adot(ss) * dy.get<double>(prefix("alpha"))
      + R_->d_ratep_d_g(ss).contract(dg.get<Symmetric>(prefix("alpha")));
  // R
  res.get<double>(dprefix("R","R")) =
        R_->d_ratep_d_h(ss)
      + R_->d_ratep_d_adot(ss) * dy.get<double>(prefix("R"))
      + R_->d_ratep_d_g(ss).contract(dg.get<Symmetric>(prefix("R")));
  // D
  res.get<double>(dprefix("R", "D")) =
        R_->d_ratep_d_D(ss)
      + R_->d_ratep_d_adot(ss) * dy.get<double>(prefix("D"))
      + R_->d_ratep_d_g(ss).contract(dg.get<Symmetric>(prefix("D")));
  // backstresses
  for (auto X: X_) {
    res.get<Symmetric>(dprefix("R", X->name())) = 
          R_->d_ratep_d_adot(ss) * dy.get<Symmetric>(prefix(X->name()))
        + dg.get<SymSymR4>(prefix(X->name())).dot(R_->d_ratep_d_g(ss)).transpose();
  }

  // D
  ss.h = state.h.get<double>(prefix("D"));
  // a
  res.get<double>(dprefix("D","alpha")) =
        D_->d_ratep_d_a(ss)
      + D_->d_ratep_d_adot(ss) * dy.get<double>(prefix("alpha"))
      + D_->d_ratep_d_g(ss).contract(dg.get<Symmetric>(prefix("alpha")));
  // R
  res.get<double>(dprefix("D","R")) =
        D_->d_ratep_d_adot(ss) * dy.get<double>(prefix("R"))
      + D_->d_ratep_d_g(ss).contract(dg.get<Symmetric>(prefix("R")));
  // D
  res.get<double>(dprefix("D","D")) =
        D_->d_ratep_d_h(ss)
      + D_->d_ratep_d_adot(ss) * dy.get<double>(prefix("D"))
      + D_->d_ratep_d_g(ss).contract(dg.get<Symmetric>(prefix("D")));
  // backstresses
  for (auto X: X_) {
    res.get<Symmetric>(dprefix("D", X->name())) = 
          D_->d_ratep_d_adot(ss) * dy.get<Symmetric>(prefix(X->name()))
        + dg.get<SymSymR4>(X->name()).dot(D_->d_ratep_d_g(ss)).transpose();
  }

  // And the backstresses...
  auto Ss = symmetric_state_(state);
  for (auto X : X_) {
    Ss.h.copy_data(state.h.get<Symmetric>(prefix(X->name())).data()); 
    // a
    res.get<Symmetric>(dprefix(X->name(),"alpha")) = 
          X->d_ratep_d_a(Ss)
        + X->d_ratep_d_adot(Ss) * dy.get<double>(prefix("alpha"))
        + X->d_ratep_d_g(Ss).dot(dg.get<Symmetric>(prefix("alpha")));

    // R
    res.get<Symmetric>(dprefix(X->name(), "R")) = 
        X->d_ratep_d_adot(Ss) * dy.get<double>(prefix("R"))
        + X->d_ratep_d_g(Ss).dot(dg.get<Symmetric>(prefix("R")));

    // D
    res.get<Symmetric>(dprefix(X->name(), "D")) = 
          X->d_ratep_d_D(Ss)
        + X->d_ratep_d_adot(Ss) * dy.get<double>(prefix("D"))
        + X->d_ratep_d_g(Ss).dot(dg.get<Symmetric>(prefix("D")));

    // backstresses
    for (auto Y : X_) {
      if (prefix(Y->name()) == prefix(X->name())) {
        res.get<SymSymR4>(dprefix(X->name(), Y->name())) = X->d_ratep_d_h(Ss);
      }
      res.get<SymSymR4>(dprefix(X->name() ,Y->name())) += 
            douter(X->d_ratep_d_adot(Ss), dy.get<Symmetric>(prefix(Y->name())))
          + X->d_ratep_d_g(Ss).dot(dg.get<SymSymR4>(prefix(Y->name())));
    }
  }
}

void WalkerFlowRule::h_time(const State & state, History & res) const
{
  // Scalar variables
  res.get<double>(prefix("alpha")) = 0.0;
  auto ss = scalar_state_(state);
  ss.h = state.h.get<double>(prefix("R"));
  res.get<double>(prefix("R")) = R_->ratet(ss);
  ss.h = state.h.get<double>(prefix("D"));
  res.get<double>(prefix("D")) = D_->ratet(ss);
  
  // Backstresses
  auto Ss = symmetric_state_(state);
  for (auto X : X_) {
    Ss.h.copy_data(state.h.get<Symmetric>(prefix(X->name())).data());
    res.get<Symmetric>(prefix(X->name())) = X->ratet(Ss);
  }
}

void WalkerFlowRule::dh_ds_time(const State & state, History & res) const
{
  // Again highly annoying, but at least there's a pattern here
  // Alpha
  res.get<Symmetric>(prefix("alpha")) = Symmetric::zero();
  
  // Common junk
  Symmetric dy;
  dy_ds(state, dy);
  SymSymR4 dg;
  dg_ds(state, dg);

  auto ss = scalar_state_(state);  

  // R
  ss.h = state.h.get<double>(prefix("R"));
  res.get<Symmetric>(prefix("R")) = 
        R_->d_ratet_d_s(ss) 
      + R_->d_ratet_d_adot(ss) * dy
      + dg.dot(R_->d_ratet_d_g(ss)).transpose();

  // D
  ss.h = state.h.get<double>(prefix("D"));
  res.get<Symmetric>(prefix("D")) = 
        D_->d_ratet_d_s(ss) 
      + D_->d_ratet_d_adot(ss) * dy
      + dg.dot(D_->d_ratet_d_g(ss)).transpose();

  // Backstresses
  auto Ss = symmetric_state_(state);
  for (auto X : X_) {
    Ss.h.copy_data(state.h.get<Symmetric>(prefix(X->name())).data());
    res.get<SymSymR4>(prefix(X->name())) = 
          X->d_ratet_d_s(Ss)
        + douter(X->d_ratet_d_adot(Ss), dy)
        + X->d_ratet_d_g(Ss).dot(dg);
  }
}

void WalkerFlowRule::dh_da_time(const State & state, History & res) const
{
  // Well this is the worst...
  
  // Common stuff we need...
  History dy = blank_derivative_<double>();
  dy_da(state, dy);

  History dg = blank_derivative_<Symmetric>();
  dg_da(state, dg);

  // Make sure the result starts out as zero
  res.zero();

  // Scalar variables
  auto ss = scalar_state_(state);

  // Alpha
  // (zero)

  // R
  ss.h = state.h.get<double>(prefix("R"));
  // a
  res.get<double>(dprefix("R","alpha")) =
        R_->d_ratet_d_a(ss)
      + R_->d_ratet_d_adot(ss) * dy.get<double>(prefix("alpha"))
      + R_->d_ratet_d_g(ss).contract(dg.get<Symmetric>(prefix("alpha")));
  // R
  res.get<double>(dprefix("R","R")) =
        R_->d_ratet_d_h(ss)
      + R_->d_ratet_d_adot(ss) * dy.get<double>(prefix("R"))
      + R_->d_ratet_d_g(ss).contract(dg.get<Symmetric>(prefix("R")));
  // D
  res.get<double>(dprefix("R", "D")) =
        R_->d_ratet_d_D(ss)
      + R_->d_ratet_d_adot(ss) * dy.get<double>(prefix("D"))
      + R_->d_ratet_d_g(ss).contract(dg.get<Symmetric>(prefix("D")));
  // backstresses
  for (auto X: X_) {
    res.get<Symmetric>(dprefix("R", X->name())) = 
          R_->d_ratet_d_adot(ss) * dy.get<Symmetric>(prefix(X->name()))
        + dg.get<SymSymR4>(prefix(X->name())).dot(R_->d_ratet_d_g(ss)).transpose();
  }

  // D
  ss.h = state.h.get<double>(prefix("D"));
  // a
  res.get<double>(dprefix("D","alpha")) =
        D_->d_ratet_d_a(ss)
      + D_->d_ratet_d_adot(ss) * dy.get<double>(prefix("alpha"))
      + D_->d_ratet_d_g(ss).contract(dg.get<Symmetric>(prefix("alpha")));
  // R
  res.get<double>(dprefix("D","R")) =
        D_->d_ratet_d_adot(ss) * dy.get<double>(prefix("R"))
      + D_->d_ratet_d_g(ss).contract(dg.get<Symmetric>(prefix("R")));
  // D
  res.get<double>(dprefix("D","D")) =
        D_->d_ratet_d_h(ss)
      + D_->d_ratet_d_adot(ss) * dy.get<double>(prefix("D"))
      + D_->d_ratet_d_g(ss).contract(dg.get<Symmetric>(prefix("D")));
  // backstresses
  for (auto X: X_) {
    res.get<Symmetric>(dprefix("D",X->name())) = 
          D_->d_ratet_d_adot(ss) * dy.get<Symmetric>(prefix(X->name()))
        + dg.get<SymSymR4>(prefix(X->name())).dot(D_->d_ratet_d_g(ss)).transpose();
  }

  // And the backstresses...
  auto Ss = symmetric_state_(state);
  for (auto X : X_) {
    Ss.h.copy_data(state.h.get<Symmetric>(prefix(X->name())).data()); 
    // a
    res.get<Symmetric>(dprefix(X->name(), "alpha")) = 
          X->d_ratet_d_a(Ss)
        + X->d_ratet_d_adot(Ss) * dy.get<double>(prefix("alpha"))
        + X->d_ratet_d_g(Ss).dot(dg.get<Symmetric>(prefix("alpha")));

    // R
    res.get<Symmetric>(dprefix(X->name(), "R")) = 
        X->d_ratet_d_adot(Ss) * dy.get<double>(prefix("R"))
        + X->d_ratet_d_g(Ss).dot(dg.get<Symmetric>(prefix("R")));

    // D
    res.get<Symmetric>(dprefix(X->name(), "D")) = 
          X->d_ratet_d_D(Ss)
        + X->d_ratet_d_adot(Ss) * dy.get<double>(prefix("D"))
        + X->d_ratet_d_g(Ss).dot(dg.get<Symmetric>(prefix("D")));

    // backstresses
    for (auto Y : X_) {
      if (prefix(Y->name()) == prefix(X->name())) {
        res.get<SymSymR4>(dprefix(X->name(), Y->name())) = X->d_ratet_d_h(Ss);
      }
      res.get<SymSymR4>(dprefix(X->name(), Y->name())) += 
            douter(X->d_ratet_d_adot(Ss), dy.get<Symmetric>(prefix(Y->name())))
          + X->d_ratet_d_g(Ss).dot(dg.get<SymSymR4>(prefix(Y->name())));
    }
  }
}

Symmetric WalkerFlowRule::TX_(const State & state) const
{
  Symmetric res = Symmetric::zero();
  for (auto X : X_)
    res += state.h.get<Symmetric>(prefix(X->name()));
  return res;
}

ScalarInternalVariable::VariableState WalkerFlowRule::scalar_state_(
    const State & state) const
{
  ScalarInternalVariable::VariableState vstate;

  vstate.a = state.h.get<double>(prefix("alpha"));
  y(state, vstate.adot);
  vstate.D = state.h.get<double>(prefix("D"));
  vstate.s = state.S;
  g(state, vstate.g);
  vstate.T = state.T;

  return vstate;
}

SymmetricInternalVariable::VariableState WalkerFlowRule::symmetric_state_(
    const State & state) const
{
  SymmetricInternalVariable::VariableState vstate;
  
  vstate.a = state.h.get<double>(prefix("alpha"));
  y(state, vstate.adot);
  vstate.D = state.h.get<double>(prefix("D"));
  vstate.s = state.S;
  g(state, vstate.g);
  vstate.T = state.T;

  return vstate;
}

double WalkerFlowRule::prefactor_(const State & state) const
{
  return eps0_->value(state.T) * 
      softening_->phi(state.h.get<double>(prefix("alpha")), state.T) * 
      scaling_->value(state.T);
}

double WalkerFlowRule::flow_(const State & state) const
{
  Symmetric d = state.S.dev() - TX_(state);
  double Y = Y_(state);
  double h = (std::sqrt(3.0/2.0) * d.norm() - Y) / state.h.get<double>(prefix("D"));
  if (h <= 0.0)
    return 0.0;
  else
    return std::pow(std::fabs(h), n_->value(state.T));
}

double WalkerFlowRule::dflow_(const State & state) const
{
  Symmetric d = state.S.dev() - TX_(state);
  double Y = Y_(state);
  double h = (std::sqrt(3.0/2.0) * d.norm() - Y) / state.h.get<double>(prefix("D"));
  if (h <= 0.0)
    return 0.0;
  else
    return n_->value(state.T) * std::pow(std::fabs(h), n_->value(state.T)-1.0);
}

double WalkerFlowRule::Y_(const State & state) const
{
  double xi = (state.h.get<double>(prefix("D")) - D_->D_0(state.T)) / D_->D_xi(state.T);
  if (xi < 0.0) xi = 0.0;
  return (k_->value(state.T) + state.h.get<double>(prefix("R"))) * std::pow(xi,
                                                                    m_->value(state.T));
}

SymSymR4 WalkerFlowRule::G_(const State & state) const
{
  Symmetric d = state.S.dev() - TX_(state);
  double nd = d.norm();
  if (nd == 0.0)
    return SymSymR4::id();
  else
    return std::sqrt(3.0/2.0) / nd * (SymSymR4::id() - douter(d/nd, d/nd));
}

void WalkerFlowRule::override_guess(double * const x)
{
  double sn = norm2_vec(x, 6);
  if (sn == 0.0) {
    x[0] += 1.0e-3;
  }
  if (x[6] < 1.0e-2) {
    x[6] += 1.0e-3;
  }
  return;
}

std::shared_ptr<ThermalScaling> default_scaling()
{
  ParameterSet params = ThermalScaling::parameters();

  return std::make_shared<ThermalScaling>(params);
}

}
