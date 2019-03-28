#include "general_flow.h"

#include "nemlmath.h"


#include <algorithm>

#include <iostream>

namespace neml {

int GeneralFlowRule::work_rate(const double * const s,
                                            const double * const alpha,
                                            const double * const edot, double T,
                                            double Tdot, double & p_dot)
{
  // By default don't calculate plastic work
  p_dot = 0.0;
  return 0;
}

int GeneralFlowRule::set_elastic_model(std::shared_ptr<LinearElasticModel>
                                       emodel)
{
  return 0;
}

TVPFlowRule::TVPFlowRule(std::shared_ptr<LinearElasticModel> elastic,
            std::shared_ptr<ViscoPlasticFlowRule> flow) :
    elastic_(elastic), flow_(flow)
{

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
  return neml::make_unique<TVPFlowRule>(
      params.get_object_parameter<LinearElasticModel>("elastic"),
      params.get_object_parameter<ViscoPlasticFlowRule>("flow")
      ); 
}


size_t TVPFlowRule::nhist() const
{
  return flow_->nhist();
}

int TVPFlowRule::init_hist(double * const h)
{
  return flow_->init_hist(h);
}

int TVPFlowRule::s(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const sdot)
{
  double erate[6];
  std::copy(edot, edot+6, erate);

  double temp[6];
  double yv;
  int ier = flow_->g(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;
  
  for (int i=0; i<6; i++) {
    erate[i] -= yv * temp[i];
  }

  ier = flow_->g_temp(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<6; i++) {
    erate[i] -= Tdot * temp[i];
  }

  ier = flow_->g_time(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<6; i++) {
    erate[i] -= temp[i];
  }
  
  double C[36];
  elastic_->C(T, C);

  mat_vec(C, 6, erate, 6, sdot);

  return 0;

}

int TVPFlowRule::ds_ds(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  double yv;
  int ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;

  double work[36];
  ier = flow_->dg_ds(s, alpha, T, work);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<36; i++) {
    work[i] *= -yv;
  }

  double t1[6];
  ier = flow_->g(s, alpha, T, t1);
  if (ier != SUCCESS) return ier;
  double t2[6];
  ier = flow_->dy_ds(s, alpha, T, t2);
  if (ier != SUCCESS) return ier;
  outer_update_minus(t1, 6, t2, 6, work);
  
  double t3[36];
  ier = flow_->dg_ds_temp(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<36; i++) {
    work[i] -= t3[i] * Tdot;
  }

  ier = flow_->dg_ds_time(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<36; i++) {
    work[i] -= t3[i];
  }

  elastic_->C(T, t3);

  mat_mat(6,6,6, t3, work, d_sdot);

  return 0;

}

int TVPFlowRule::ds_da(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  double yv;
  int ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;
  
  int sz = 6 * nhist();
  
  std::vector<double> workv(sz);
  double * work = &workv[0];
  ier = flow_->dg_da(s, alpha, T, work);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) {
    work[i] *= -yv;
  }

  double t1[6];
  ier = flow_->g(s, alpha, T, t1);
  if (ier != SUCCESS) return ier;
  std::vector<double> t2v(nhist());
  double * t2 = &t2v[0];
  ier = flow_->dy_da(s, alpha, T, t2);
  if (ier != SUCCESS) return ier;
  outer_update_minus(t1, 6, t2, nhist(), work);
  
  std::vector<double> t3v(sz);
  double * t3 = &t3v[0];
  ier = flow_->dg_da_temp(s, alpha, T, t3);
  if (ier != SUCCESS) return ier; 
  for (int i=0; i<sz; i++) {
    work[i] -= t3[i] * Tdot;
  }

  ier = flow_->dg_da_time(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) {
    work[i] -= t3[i];
  }

  double C[36];
  elastic_->C(T, C);

  mat_mat(6, nhist(), 6, C, work, d_sdot);

  return 0;

}

int TVPFlowRule::ds_de(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_sdot)
{
  return elastic_->C(T, d_sdot);
}

int TVPFlowRule::a(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const adot)
{
  double dg;
  int ier = flow_->y(s, alpha, T, dg);
  if (ier != SUCCESS) return 0;

  ier = flow_->h(s, alpha, T, adot);
  if (ier != SUCCESS) return 0;
  for (size_t i=0; i<nhist(); i++) adot[i] *= dg;
  
  std::vector<double> tempv(nhist());
  double * temp = &tempv[0];
  ier = flow_->h_temp(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (size_t i=0; i<nhist(); i++) adot[i] += temp[i] * Tdot;

  ier = flow_->h_time(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (size_t i=0; i<nhist(); i++) adot[i] += temp[i];

  return 0;

}

int TVPFlowRule::da_ds(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  double dg;
  int ier = flow_->y(s, alpha, T, dg);
  if (ier != SUCCESS) return ier;

  int sz = nhist() * 6;

  ier = flow_->dh_ds(s, alpha, T, d_adot);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] *= dg;

  std::vector<double> t1v(nhist());
  double * t1 = &t1v[0];
  ier = flow_->h(s, alpha, T, t1);
  if (ier != SUCCESS) return ier;

  double t2[6];
  ier = flow_->dy_ds(s, alpha, T, t2);
  if (ier != SUCCESS) return ier;

  outer_update(t1, nhist(), t2, 6, d_adot);
  
  std::vector<double> t3v(sz);
  double * t3 = &t3v[0];
  ier = flow_->dh_ds_temp(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] += t3[i] * Tdot;

  ier = flow_->dh_ds_time(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] += t3[i];

  return 0;
  
}

int TVPFlowRule::da_da(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  double dg;
  int ier = flow_->y(s, alpha, T, dg);
  if (ier != SUCCESS) return ier;

  int sz = nhist() * nhist();

  ier = flow_->dh_da(s, alpha, T, d_adot);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] *= dg;
  
  std::vector<double> t1v(nhist());
  double * t1 = &t1v[0];
  ier = flow_->h(s, alpha, T, t1);
  if (ier != SUCCESS) return ier;
  
  std::vector<double> t2v(nhist());
  double * t2 = &t2v[0];
  ier = flow_->dy_da(s, alpha, T, t2);
  if (ier != SUCCESS) return ier;

  outer_update(t1, nhist(), t2, nhist(), d_adot);
  
  std::vector<double> t3v(sz);
  double * t3 = &t3v[0];
  ier = flow_->dh_da_temp(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] += t3[i] * Tdot;

  ier = flow_->dh_da_time(s, alpha, T, t3);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<sz; i++) d_adot[i] += t3[i];

  return 0;
}

int TVPFlowRule::da_de(const double * const s, const double * const alpha,
              const double * const edot, double T,
              double Tdot,
              double * const d_adot)
{
  std::fill(d_adot, d_adot+(nhist()*6), 0.0);

  return 0;
}

int TVPFlowRule::work_rate(const double * const s,
                                    const double * const alpha,
                                    const double * const edot, double T,
                                    double Tdot, double & p_dot)
{
  double erate[6];
  std::fill(erate, erate+6, 0.0);

  double temp[6];
  double yv;
  int ier = flow_->g(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  ier = flow_->y(s, alpha, T, yv);
  if (ier != SUCCESS) return ier;

  for (int i=0; i<6; i++) {
    erate[i] += yv * temp[i];
  }

  ier = flow_->g_temp(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<6; i++) {
    erate[i] += Tdot * temp[i];
  }

  ier = flow_->g_time(s, alpha, T, temp);
  if (ier != SUCCESS) return ier;
  for (int i=0; i<6; i++) {
    erate[i] += temp[i];
  }

  p_dot = dot_vec(s, erate, 6);
  
  return 0;
}

int TVPFlowRule::elastic_strains(const double * const s_np1, double T_np1,
                                 double * const e_np1) const
{
  double S[36];
  elastic_->S(T_np1, S);
  mat_vec(S, 6, s_np1, 6, e_np1);

  return 0;
}

int TVPFlowRule::set_elastic_model(std::shared_ptr<LinearElasticModel> emodel)
{
  elastic_ = emodel;
  return 0;
}

} // namespace neml
