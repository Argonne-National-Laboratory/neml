#include "general_flow.h"

#include "math/nemlmath.h"


#include <algorithm>

#include <iostream>

namespace neml {

GeneralFlowRule::GeneralFlowRule(ParameterSet & params) :
    HistoryNEMLObject(params)
{
}

double GeneralFlowRule::work_rate(const Symmetric & s, const History & alpha, 
                                  const Symmetric & edot, double T, double Tdot)
{
  // By default don't calculate plastic work
  return 0.0;
}

void GeneralFlowRule::set_elastic_model(std::shared_ptr<LinearElasticModel>
                                       emodel)
{
}

void GeneralFlowRule::override_guess(double * const x)
{
  return;
}

TVPFlowRule::TVPFlowRule(ParameterSet & params) :
    GeneralFlowRule(params),
    elastic_(params.get_object_parameter<LinearElasticModel>("elastic")),
    flow_(params.get_object_parameter<ViscoPlasticFlowRule>("flow"))
{
  cache_history_();
}

std::string TVPFlowRule::type()
{
  return "TVPFlowRule";
}

ParameterSet TVPFlowRule::parameters()
{
  ParameterSet pset(TVPFlowRule::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("flow");

  return pset;
}

std::unique_ptr<NEMLObject> TVPFlowRule::initialize(ParameterSet & params)
{
  return neml::make_unique<TVPFlowRule>(params); 
}


void TVPFlowRule::populate_hist(History & h) const
{
  flow_->set_variable_prefix(get_variable_prefix());
  flow_->populate_hist(h);
}

void TVPFlowRule::init_hist(History & h) const
{
  flow_->init_hist(h);
}

Symmetric TVPFlowRule::s(const Symmetric & s, const History & alpha, 
                         const Symmetric & edot, double T, double Tdot)
{
  Symmetric temp;
  double yv;
  flow_->g(s.data(), alpha.rawptr(), T, temp.s());
  flow_->y(s.data(), alpha.rawptr(), T, yv);
  if (yv > NEML_STRAIN_RATE_LIMIT) 
    throw NEMLError("Strain rate exceeds rate limit");
  
  Symmetric erate = edot - yv * temp;

  flow_->g_temp(s.data(), alpha.rawptr(), T, temp.s());
  erate -= Tdot * temp;

  flow_->g_time(s.data(), alpha.rawptr(), T, temp.s());
  erate -= temp;

  return elastic_->C(T).dot(erate);
}

SymSymR4 TVPFlowRule::ds_ds(const Symmetric & s, const History & alpha, 
                            const Symmetric & edot, double T, double Tdot)
{
  double yv;
  flow_->y(s.data(), alpha.rawptr(), T, yv);

  SymSymR4 work;
  flow_->dg_ds(s.data(), alpha.rawptr(), T, work.s());
  work *= -yv;

  Symmetric t1, t2;
  flow_->g(s.data(), alpha.rawptr(), T, t1.s());
  flow_->dy_ds(s.data(), alpha.rawptr(), T, t2.s());
  work -= douter(t1, t2);

  SymSymR4 t3;
  flow_->dg_ds_temp(s.data(), alpha.rawptr(), T, t3.s());
  work -= t3 * Tdot;
  
  flow_->dg_ds_time(s.data(), alpha.rawptr(), T, t3.s());
  work -= t3;

  return elastic_->C(T).dot(work);
}

History TVPFlowRule::ds_da(const Symmetric & s, const History & alpha, 
                           const Symmetric & edot, double T, double Tdot)
{
  double yv;
  flow_->y(s.data(), alpha.rawptr(), T, yv);
  
  History work = gather_blank_history_().derivative<Symmetric>();
  flow_->dg_da(s.data(), alpha.rawptr(), T, work.rawptr());
  work *= -yv;
  
  // Ugh, we need to make these matrices...
  Symmetric t1;
  flow_->g(s.data(), alpha.rawptr(), T, t1.s());
  History t2 = gather_blank_history_();
  flow_->dy_da(s.data(), alpha.rawptr(), T, t2.rawptr());
  outer_update_minus(t1.data(), 6, t2.rawptr(), nhist(), work.rawptr());
  
  History t3 = gather_blank_history_().derivative<Symmetric>();
  flow_->dg_da_temp(s.data(), alpha.rawptr(), T, t3.rawptr());
  work -= t3 * Tdot;

  flow_->dg_da_time(s.data(), alpha.rawptr(), T, t3.rawptr());
  work -= t3;
  
  SymSymR4 C = elastic_->C(T);

  return work.premultiply(C);
}

SymSymR4 TVPFlowRule::ds_de(const Symmetric & s, const History & alpha, 
                           const Symmetric & edot, double T, double Tdot)
{
  return elastic_->C(T);
}

History TVPFlowRule::a(const Symmetric & s, const History & alpha, 
                       const Symmetric & edot, double T, double Tdot)
{
  double dg;
  flow_->y(s.data(), alpha.rawptr(), T, dg);
  
  History adot = gather_blank_history_();
  flow_->h(s.data(), alpha.rawptr(), T, adot.rawptr());
  adot *= dg;
  
  History temp = gather_blank_history_();
  flow_->h_temp(s.data(), alpha.rawptr(), T, temp.rawptr());
  adot += temp * Tdot;

  flow_->h_time(s.data(), alpha.rawptr(), T, temp.rawptr());
  adot += temp;

  return adot;
}

History TVPFlowRule::da_ds(const Symmetric & s, const History & alpha, 
                           const Symmetric & edot, double T, double Tdot)
{
  double dg;
  flow_->y(s.data(), alpha.rawptr(), T, dg);

  History d_adot = gather_blank_history_().derivative<Symmetric>();
  flow_->dh_ds(s.data(), alpha.rawptr(), T, d_adot.rawptr());
  d_adot *= dg;
  
  // Ugh, history should be a matrix
  History t1 = gather_blank_history_();
  flow_->h(s.data(), alpha.rawptr(), T, t1.rawptr());
  Symmetric t2;
  flow_->dy_ds(s.data(), alpha.rawptr(), T, t2.s());
  outer_update(t1.rawptr(), nhist(), t2.data(), 6, d_adot.rawptr());
  
  History t3 = gather_blank_history_().derivative<Symmetric>();
  flow_->dh_ds_temp(s.data(), alpha.rawptr(), T, t3.rawptr());
  d_adot += t3 * Tdot;

  flow_->dh_ds_time(s.data(), alpha.rawptr(), T, t3.rawptr());
  return d_adot + t3;
}

History TVPFlowRule::da_da(const Symmetric & s, const History & alpha, 
                           const Symmetric & edot, double T, double Tdot)
{
  double dg;
  flow_->y(s.data(), alpha.rawptr(), T, dg);

  int nh = nhist();
  
  History d_adot = gather_blank_history_().derivative<History>();
  flow_->dh_da(s.data(), alpha.rawptr(), T, d_adot.rawptr());
  d_adot *= dg;
  
  // Ugh
  History t1 = gather_blank_history_();
  flow_->h(s.data(), alpha.rawptr(), T, t1.rawptr());
  History t2 = gather_blank_history_();
  flow_->dy_da(s.data(), alpha.rawptr(), T, t2.rawptr());
  outer_update(t1.rawptr(), nh, t2.rawptr(), nh, d_adot.rawptr());
  
  History t3 = gather_blank_history_().derivative<History>();
  flow_->dh_da_temp(s.data(), alpha.rawptr(), T, t3.rawptr());
  d_adot += t3 * Tdot;
  
  flow_->dh_da_time(s.data(), alpha.rawptr(), T, t3.rawptr());
  return d_adot + t3;
}

History TVPFlowRule::da_de(const Symmetric & s, const History & alpha, 
                           const Symmetric & edot, double T, double Tdot)
{
  return gather_blank_history_().derivative<Symmetric>();
}

double TVPFlowRule::work_rate(const Symmetric & s, const History & alpha, 
                              const Symmetric & edot, double T, double Tdot)
{
  Symmetric dir;
  double yv;
  flow_->g(s.data(), alpha.rawptr(), T, dir.s());
  flow_->y(s.data(), alpha.rawptr(), T, yv);

  Symmetric erate = yv * dir;

  flow_->g_temp(s.data(), alpha.rawptr(), T, dir.s());
  erate += Tdot * dir;

  flow_->g_time(s.data(), alpha.rawptr(), T, dir.s());
  erate += dir;
  
  return s.contract(erate);
}

void TVPFlowRule::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
}

void TVPFlowRule::override_guess(double * const x)
{
  flow_->override_guess(x);
}

} // namespace neml
