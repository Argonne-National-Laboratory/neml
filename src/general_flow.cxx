#include "general_flow.h"

#include "math/nemlmath.h"


#include <algorithm>

#include <iostream>

namespace neml {

GeneralFlowRule::GeneralFlowRule(ParameterSet & params) :
    HistoryNEMLObject(params)
{
}

void GeneralFlowRule::work_rate(const double * const s,
                                            const double * const alpha,
                                            const double * const edot, double T,
                                            double Tdot, double & p_dot)
{
  // By default don't calculate plastic work
  p_dot = 0.0;
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

void TVPFlowRule::s(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const sdot)
{
  double erate[6];
  std::copy(edot, edot+6, erate);

  double temp[6];
  double yv;
  flow_->g(s, alpha, T, temp);
  flow_->y(s, alpha, T, yv);
  if (yv > NEML_STRAIN_RATE_LIMIT) 
    throw NEMLError("Strain rate exceeds rate limit");
  
  for (int i=0; i<6; i++) {
    erate[i] -= yv * temp[i];
  }

  flow_->g_temp(s, alpha, T, temp);
  for (int i=0; i<6; i++) {
    erate[i] -= Tdot * temp[i];
  }

  flow_->g_time(s, alpha, T, temp);
  for (int i=0; i<6; i++) {
    erate[i] -= temp[i];
  }
  
  double C[36];
  elastic_->C(T, C);

  mat_vec(C, 6, erate, 6, sdot);


}

void TVPFlowRule::ds_ds(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  double yv;
  flow_->y(s, alpha, T, yv);

  double work[36];
  flow_->dg_ds(s, alpha, T, work);
  for (int i=0; i<36; i++) {
    work[i] *= -yv;
  }

  double t1[6];
  flow_->g(s, alpha, T, t1);
  double t2[6];
  flow_->dy_ds(s, alpha, T, t2);
  outer_update_minus(t1, 6, t2, 6, work);
  
  double t3[36];
  flow_->dg_ds_temp(s, alpha, T, t3);
  for (int i=0; i<36; i++) {
    work[i] -= t3[i] * Tdot;
  }

  flow_->dg_ds_time(s, alpha, T, t3);
  for (int i=0; i<36; i++) {
    work[i] -= t3[i];
  }

  elastic_->C(T, t3);

  mat_mat(6,6,6, t3, work, d_sdot);


}

void TVPFlowRule::ds_da(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  double yv;
  flow_->y(s, alpha, T, yv);
  
  int sz = 6 * nh();
  
  std::vector<double> workv(sz);
  double * work = &workv[0];
  flow_->dg_da(s, alpha, T, work);
  for (int i=0; i<sz; i++) {
    work[i] *= -yv;
  }

  double t1[6];
  flow_->g(s, alpha, T, t1);
  std::vector<double> t2v(nh());
  double * t2 = &t2v[0];
  flow_->dy_da(s, alpha, T, t2);
  outer_update_minus(t1, 6, t2, nh(), work);
  
  std::vector<double> t3v(sz);
  double * t3 = &t3v[0];
  flow_->dg_da_temp(s, alpha, T, t3);
  for (int i=0; i<sz; i++) {
    work[i] -= t3[i] * Tdot;
  }

  flow_->dg_da_time(s, alpha, T, t3);
  for (int i=0; i<sz; i++) {
    work[i] -= t3[i];
  }

  double C[36];
  elastic_->C(T, C);

  mat_mat(6, nh(), 6, C, work, d_sdot);
}

void TVPFlowRule::ds_de(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  elastic_->C(T, d_sdot);
}

void TVPFlowRule::a(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const adot)
{
  double dg;
  flow_->y(s, alpha, T, dg);

  flow_->h(s, alpha, T, adot);
  for (size_t i=0; i<nh(); i++) adot[i] *= dg;
  
  std::vector<double> tempv(nh());
  double * temp = &tempv[0];
  flow_->h_temp(s, alpha, T, temp);
  for (size_t i=0; i<nh(); i++) adot[i] += temp[i] * Tdot;

  flow_->h_time(s, alpha, T, temp);
  for (size_t i=0; i<nh(); i++) adot[i] += temp[i];


}

void TVPFlowRule::da_ds(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  double dg;
  flow_->y(s, alpha, T, dg);

  int sz = nh() * 6;

  flow_->dh_ds(s, alpha, T, d_adot);
  for (int i=0; i<sz; i++) d_adot[i] *= dg;

  std::vector<double> t1v(nh());
  double * t1 = &t1v[0];
  flow_->h(s, alpha, T, t1);

  double t2[6];
  flow_->dy_ds(s, alpha, T, t2);

  outer_update(t1, nh(), t2, 6, d_adot);
  
  std::vector<double> t3v(sz);
  double * t3 = &t3v[0];
  flow_->dh_ds_temp(s, alpha, T, t3);
  for (int i=0; i<sz; i++) d_adot[i] += t3[i] * Tdot;

  flow_->dh_ds_time(s, alpha, T, t3);
  for (int i=0; i<sz; i++) d_adot[i] += t3[i];

  
}

void TVPFlowRule::da_da(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  double dg;
  flow_->y(s, alpha, T, dg);

  int nhi = nh();
  int sz = nhi * nhi;

  flow_->dh_da(s, alpha, T, d_adot);
  for (int i=0; i<sz; i++) d_adot[i] *= dg;
  
  std::vector<double> t1v(nhi);
  double * t1 = &t1v[0];
  flow_->h(s, alpha, T, t1);
  
  std::vector<double> t2v(nhi);
  double * t2 = &t2v[0];
  flow_->dy_da(s, alpha, T, t2);

  outer_update(t1, nhi, t2, nhi, d_adot);
  
  std::vector<double> t3v(sz);
  double * t3 = &t3v[0];
  flow_->dh_da_temp(s, alpha, T, t3);
  for (int i=0; i<sz; i++) d_adot[i] += t3[i] * Tdot;

  flow_->dh_da_time(s, alpha, T, t3);
  for (int i=0; i<sz; i++) d_adot[i] += t3[i];

}

void TVPFlowRule::da_de(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  std::fill(d_adot, d_adot+(nh()*6), 0.0);

}

void TVPFlowRule::work_rate(const double * const s,
                                    const double * const alpha,
                                    const double * const edot, double T,
                                    double Tdot, double & p_dot)
{
  double erate[6];
  std::fill(erate, erate+6, 0.0);

  double temp[6];
  double yv;
  flow_->g(s, alpha, T, temp);
  flow_->y(s, alpha, T, yv);

  for (int i=0; i<6; i++) {
    erate[i] += yv * temp[i];
  }

  flow_->g_temp(s, alpha, T, temp);
  for (int i=0; i<6; i++) {
    erate[i] += Tdot * temp[i];
  }

  flow_->g_time(s, alpha, T, temp);
  for (int i=0; i<6; i++) {
    erate[i] += temp[i];
  }

  p_dot = dot_vec(s, erate, 6);
}

void TVPFlowRule::elastic_strains(const double * const s_np1, double T_np1,
                                 double * const e_np1) const
{
  double S[36];
  elastic_->S(T_np1, S);
  mat_vec(S, 6, s_np1, 6, e_np1);

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
