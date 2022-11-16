#include "visco_flow.h"

#include "math/nemlmath.h"

#include <cmath>
#include <iostream>

namespace neml {

ViscoPlasticFlowRule::ViscoPlasticFlowRule(ParameterSet & params) :
    HistoryNEMLObject(params)
{

}

// Default implementation of flow rule wrt time
void ViscoPlasticFlowRule::g_time(const double * const s, 
                                 const double * const alpha, double T, 
                                 double * const gv) const
{
  std::fill(gv, gv+6, 0.0);
}

void ViscoPlasticFlowRule::dg_ds_time(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dgv) const
{
  std::fill(dgv, dgv+36, 0.0);
}

void ViscoPlasticFlowRule::dg_da_time(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dgv) const
{
  std::fill(dgv, dgv+6*nh(), 0.0);
}

// Default implementation of flow rule wrt temperature
void ViscoPlasticFlowRule::g_temp(const double * const s, 
                                 const double * const alpha, double T, 
                                 double * const gv) const
{
  std::fill(gv, gv+6, 0.0);
}

void ViscoPlasticFlowRule::dg_ds_temp(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dgv) const
{
  std::fill(dgv, dgv+36, 0.0);
}

void ViscoPlasticFlowRule::dg_da_temp(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dgv) const
{
  std::fill(dgv, dgv+6*nh(), 0.0);
}

// Default implementation of hardening rule wrt time
void ViscoPlasticFlowRule::h_time(const double * const s, 
                                 const double * const alpha, double T, 
                                 double * const hv) const
{
  std::fill(hv, hv+nh(), 0.0);
}

void ViscoPlasticFlowRule::dh_ds_time(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dhv) const
{
  std::fill(dhv, dhv+6*nh(), 0.0);
}

void ViscoPlasticFlowRule::dh_da_time(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dhv) const
{
  std::fill(dhv, dhv+nh()*nh(), 0.0);
}

// Default implementation of hardening rule wrt temperature
void ViscoPlasticFlowRule::h_temp(const double * const s, 
                                 const double * const alpha, double T, 
                                 double * const hv) const
{
  std::fill(hv, hv+nh(), 0.0);
}

void ViscoPlasticFlowRule::dh_ds_temp(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dhv) const
{
  std::fill(dhv, dhv+6*nh(), 0.0);
}

void ViscoPlasticFlowRule::dh_da_temp(const double * const s, 
                                     const double * const alpha, double T,
                                     double * const dhv) const
{
  std::fill(dhv, dhv+nh()*nh(), 0.0);
}

void ViscoPlasticFlowRule::override_guess(double * const guess)
{
  return;
}

SuperimposedViscoPlasticFlowRule::SuperimposedViscoPlasticFlowRule(ParameterSet & params) :
    ViscoPlasticFlowRule(params),
    rules_(params.get_object_parameter_vector<ViscoPlasticFlowRule>("flow_rules")),
    offsets_(rules_.size()+1)
{
  offsets_[0] = 0;
  for (size_t i = 1; i <= rules_.size(); i++) {
    offsets_[i] = offsets_[i-1] + rules_[i-1]->nh();
  }

  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->set_variable_prefix("model"+std::to_string(i)+"_");

  cache_history_();
}

std::string SuperimposedViscoPlasticFlowRule::type()
{
  return "SuperimposedViscoPlasticFlowRule";
}

ParameterSet SuperimposedViscoPlasticFlowRule::parameters()
{
  ParameterSet pset(SuperimposedViscoPlasticFlowRule::type());

  pset.add_parameter<std::vector<NEMLObject>>("flow_rules");

  return pset;
}

std::unique_ptr<NEMLObject> SuperimposedViscoPlasticFlowRule::initialize(
    ParameterSet & params)
{
  return neml::make_unique<SuperimposedViscoPlasticFlowRule>(params); 
}

size_t SuperimposedViscoPlasticFlowRule::nmodels() const
{
  return rules_.size();
}

void SuperimposedViscoPlasticFlowRule::populate_hist(History & hist) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->populate_hist(hist);
}

void SuperimposedViscoPlasticFlowRule::init_hist(History & hist) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->init_hist(hist);
}

void SuperimposedViscoPlasticFlowRule::y(const double* const s, 
                                         const double* const alpha, double T,
                                         double & yv) const
{
  yv = 0.0;
  double yvi;
  for (size_t i = 0; i < nmodels(); i++) {
    rules_[i]->y(s, model_history_(alpha, i), T, yvi);
    yv += yvi;
  }
}

void SuperimposedViscoPlasticFlowRule::dy_ds(const double* const s, 
                                             const double* const alpha, double T,
                                             double * const dyv) const
{
  std::fill(dyv, dyv+6, 0.0);
  double dyvi[6];
  for (size_t i = 0; i < nmodels(); i++) {
    rules_[i]->dy_ds(s, model_history_(alpha, i), T, dyvi);
    add_vec(dyv, dyvi, 6, dyv);
  }
}

void SuperimposedViscoPlasticFlowRule::dy_da(const double* const s, 
                                             const double* const alpha, double T,
                                             double * const dyv) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->dy_da(s, model_history_(alpha, i), T, model_history_(dyv, i));
}

void SuperimposedViscoPlasticFlowRule::g(const double * const s, 
                                         const double * const alpha, double T,
                                         double * const gv) const
{
  std::fill(gv, gv+6, 0.0);
  double ysum = 0.0;
  double yi;
  double gvi[6];
  for (size_t i = 0; i < nmodels(); i++) {
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g(s, model_history_(alpha, i), T, gvi);
    for (size_t j = 0; j < 6; j++) gvi[j] *= yi;
    ysum += yi;
    add_vec(gv, gvi, 6, gv);
  }
  if (ysum > 0)
    for (size_t j = 0; j < 6; j++) gv[j] /= ysum;
}

void SuperimposedViscoPlasticFlowRule::dg_ds(const double * const s, 
                                             const double * const alpha, 
                                             double T,
                                             double * const dgv) const
{
  double yv;
  y(s, alpha, T, yv); 
  
  double yi;
  double dgvi[36];
  double dyvi[6];
  double gi[6];

  std::fill(dgv, dgv+36, 0.0);

  for (size_t i = 0; i < nmodels(); i++) {
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g(s, model_history_(alpha, i), T, gi);
    rules_[i]->dg_ds(s, model_history_(alpha, i), T, dgvi);
    rules_[i]->dy_ds(s, model_history_(alpha, i), T, dyvi);
    for (size_t j = 0; j < 36; j++) dgv[j] += yi * dgvi[j];
    outer_update(gi, 6, dyvi, 6, dgv);
  }
  if (yv > 0)
    for (size_t i = 0; i < 36; i++) dgv[i] /= yv;
  
  double dy[6];
  dy_ds(s, alpha, T, dy);
  double gv[6];
  g(s, alpha, T, gv);
  if (yv > 0)
    for (size_t i = 0; i < 6; i++) gv[i] /= yv;

  outer_update_minus(gv, 6, dy, 6, dgv);
}

void SuperimposedViscoPlasticFlowRule::dg_da(const double * const s, 
                                             const double * const alpha, 
                                             double T,
                                             double * const dgv) const
{
  std::fill(dgv, dgv+6*nh(), 0.0);
  
  double yv;
  y(s, alpha, T, yv);
  
  double yi;
  double gi[6];

  for (size_t i = 0; i < nmodels(); i++) {
    double * dgvi = new double [6*rules_[i]->nh()];
    double * dyi = new double [rules_[i]->nh()]; 
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g(s, model_history_(alpha, i), T, gi);
    rules_[i]->dg_da(s, model_history_(alpha,i), T, dgvi);
    rules_[i]->dy_da(s, model_history_(alpha, i), T, dyi);

    for (size_t a = 0; a < 6; a++) {
      for (size_t b = 0; b < rules_[i]->nh(); b++) {
        dgv[CINDEX(a,b+offsets_[i],nh())] =
            yi*dgvi[CINDEX(a,b,rules_[i]->nh())] + 
            gi[a] * dyi[b];
      }
    }

    delete [] dgvi;
    delete [] dyi;
  }
  
  if (yv > 0)
    for (size_t i = 0; i < 6*nh(); i++) dgv[i] /= yv;

  double * dy = new double [nh()];
  double gv[6];
  dy_da(s, alpha, T, dy);
  g(s, alpha, T, gv);
  
  if (yv > 0)
    for (size_t i = 0; i < 6; i++) gv[i] /= yv;

  outer_update_minus(gv, 6, dy, nh(), dgv);

  delete [] dy;
}

void SuperimposedViscoPlasticFlowRule::g_time(const double * const s, 
                                         const double * const alpha, double T,
                                         double * const gv) const
{
  std::fill(gv, gv+6, 0.0);
  double ysum = 0.0;
  double yi;
  double gvi[6];
  for (size_t i = 0; i < nmodels(); i++) {
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g_time(s, model_history_(alpha, i), T, gvi);
    for (size_t j = 0; j < 6; j++) gvi[j] *= yi;
    ysum += yi;
    add_vec(gv, gvi, 6, gv);
  }
  if (ysum > 0)
    for (size_t j = 0; j < 6; j++) gv[j] /= ysum;
}

void SuperimposedViscoPlasticFlowRule::dg_ds_time(const double * const s, 
                                             const double * const alpha, 
                                             double T,
                                             double * const dgv) const
{
  double yv;
  y(s, alpha, T, yv); 
  
  double yi;
  double dgvi[36];
  double dyvi[6];
  double gi[6];

  std::fill(dgv, dgv+36, 0.0);

  for (size_t i = 0; i < nmodels(); i++) {
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g_time(s, model_history_(alpha, i), T, gi);
    rules_[i]->dg_ds_time(s, model_history_(alpha, i), T, dgvi);
    rules_[i]->dy_ds(s, model_history_(alpha, i), T, dyvi);
    for (size_t j = 0; j < 36; j++) dgv[j] += yi * dgvi[j];
    outer_update(gi, 6, dyvi, 6, dgv);
  }
  if (yv > 0)
    for (size_t i = 0; i < 36; i++) dgv[i] /= yv;
  
  double dy[6];
  dy_ds(s, alpha, T, dy);
  double gv[6];
  g_time(s, alpha, T, gv);

  if (yv > 0)
    for (size_t i = 0; i < 6; i++) gv[i] /= yv;

  outer_update_minus(gv, 6, dy, 6, dgv);
}

void SuperimposedViscoPlasticFlowRule::dg_da_time(const double * const s, 
                                             const double * const alpha, 
                                             double T,
                                             double * const dgv) const
{
  std::fill(dgv, dgv+6*nh(), 0.0);
  
  double yv;
  y(s, alpha, T, yv);
  
  double yi;
  double gi[6];

  for (size_t i = 0; i < nmodels(); i++) {
    double * dgvi = new double [6*rules_[i]->nh()];
    double * dyi = new double [rules_[i]->nh()]; 
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g_time(s, model_history_(alpha, i), T, gi);
    rules_[i]->dg_da_time(s, model_history_(alpha,i), T, dgvi);
    rules_[i]->dy_da(s, model_history_(alpha, i), T, dyi);

    for (size_t a = 0; a < 6; a++) {
      for (size_t b = 0; b < rules_[i]->nh(); b++) {
        dgv[CINDEX(a,b+offsets_[i],nh())] =
            yi*dgvi[CINDEX(a,b,rules_[i]->nh())] + 
            gi[a] * dyi[b];
      }
    }

    delete [] dgvi;
    delete [] dyi;
  }
  
  if (yv > 0)
    for (size_t i = 0; i < 6*nh(); i++) dgv[i] /= yv;

  double * dy = new double [nh()];
  double gv[6];
  dy_da(s, alpha, T, dy);
  g_time(s, alpha, T, gv);
  
  if (yv > 0)
    for (size_t i = 0; i < 6; i++) gv[i] /= yv;

  outer_update_minus(gv, 6, dy, nh(), dgv);

  delete [] dy;
}

void SuperimposedViscoPlasticFlowRule::g_temp(const double * const s, 
                                         const double * const alpha, double T,
                                         double * const gv) const
{
  std::fill(gv, gv+6, 0.0);
  double ysum = 0.0;
  double yi;
  double gvi[6];
  for (size_t i = 0; i < nmodels(); i++) {
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g_temp(s, model_history_(alpha, i), T, gvi);
    for (size_t j = 0; j < 6; j++) gvi[j] *= yi;
    ysum += yi;
    add_vec(gv, gvi, 6, gv);
  }
  if (ysum > 0)
    for (size_t j = 0; j < 6; j++) gv[j] /= ysum;
}

void SuperimposedViscoPlasticFlowRule::dg_ds_temp(const double * const s, 
                                             const double * const alpha, 
                                             double T,
                                             double * const dgv) const
{
  double yv;
  y(s, alpha, T, yv); 
  
  double yi;
  double dgvi[36];
  double dyvi[6];
  double gi[6];

  std::fill(dgv, dgv+36, 0.0);

  for (size_t i = 0; i < nmodels(); i++) {
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g_temp(s, model_history_(alpha, i), T, gi);
    rules_[i]->dg_ds_temp(s, model_history_(alpha, i), T, dgvi);
    rules_[i]->dy_ds(s, model_history_(alpha, i), T, dyvi);
    for (size_t j = 0; j < 36; j++) dgv[j] += yi * dgvi[j];
    outer_update(gi, 6, dyvi, 6, dgv);
  }
  if (yv > 0)
    for (size_t i = 0; i < 36; i++) dgv[i] /= yv;
  
  double dy[6];
  dy_ds(s, alpha, T, dy);
  double gv[6];
  g_temp(s, alpha, T, gv);
  if (yv > 0)
    for (size_t i = 0; i < 6; i++) gv[i] /= yv;

  outer_update_minus(gv, 6, dy, 6, dgv);
}

void SuperimposedViscoPlasticFlowRule::dg_da_temp(const double * const s, 
                                             const double * const alpha, 
                                             double T,
                                             double * const dgv) const
{
  std::fill(dgv, dgv+6*nh(), 0.0);
  
  double yv;
  y(s, alpha, T, yv);
  
  double yi;
  double gi[6];

  for (size_t i = 0; i < nmodels(); i++) {
    double * dgvi = new double [6*rules_[i]->nh()];
    double * dyi = new double [rules_[i]->nh()]; 
    rules_[i]->y(s, model_history_(alpha, i), T, yi);
    rules_[i]->g_temp(s, model_history_(alpha, i), T, gi);
    rules_[i]->dg_da_temp(s, model_history_(alpha,i), T, dgvi);
    rules_[i]->dy_da(s, model_history_(alpha, i), T, dyi);

    for (size_t a = 0; a < 6; a++) {
      for (size_t b = 0; b < rules_[i]->nh(); b++) {
        dgv[CINDEX(a,b+offsets_[i],nh())] =
            yi*dgvi[CINDEX(a,b,rules_[i]->nh())] + 
            gi[a] * dyi[b];
      }
    }

    delete [] dgvi;
    delete [] dyi;
  }
  
  if (yv > 0)
    for (size_t i = 0; i < 6*nh(); i++) dgv[i] /= yv;

  double * dy = new double [nh()];
  double gv[6];
  dy_da(s, alpha, T, dy);
  g_temp(s, alpha, T, gv);
  
  if (yv > 0)
    for (size_t i = 0; i < 6; i++) gv[i] /= yv;

  outer_update_minus(gv, 6, dy, nh(), dgv);

  delete [] dy;
}

void SuperimposedViscoPlasticFlowRule::h(const double * const s, 
                                         const double * const alpha, double T,
                                         double * const hv) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->h(s, model_history_(alpha,i), T, model_history_(hv,i));
}

void SuperimposedViscoPlasticFlowRule::dh_ds(const double * const s,
                                             const double * const alpha,
                                             double T,
                                             double * const dhv) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->dh_ds(s, model_history_(alpha,i), T, &(dhv[offsets_[i]*6]));
}

void SuperimposedViscoPlasticFlowRule::dh_da(const double * const s, 
                                             const double * const alpha, 
                                             double T,
                                             double * const dhv) const
{
  for (size_t i = 0; i < nmodels(); i++) {
    size_t nhi = rules_[i]->nh();
    double * dhvi = new double [nhi*nhi];
    rules_[i]->dh_da(s, model_history_(alpha,i), T, dhvi);
    for (size_t a = 0; a < nhi; a++)
      for (size_t b = 0; b < nhi; b++)
        dhv[CINDEX(a+offsets_[i],b+offsets_[i],nh())] = dhvi[CINDEX(a,b,nhi)];

    delete [] dhvi;
  }
}

void SuperimposedViscoPlasticFlowRule::h_time(const double * const s, 
                                              const double * const alpha,
                                              double T,
                                              double * const hv) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->h_time(s, model_history_(alpha,i), T, model_history_(hv,i));
}

void SuperimposedViscoPlasticFlowRule::dh_ds_time(const double * const s, 
                                                  const double * const alpha, 
                                                  double T,
                                                  double * const dhv) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->dh_ds_time(s, model_history_(alpha,i), T, &(dhv[offsets_[i]*6]));
}

void SuperimposedViscoPlasticFlowRule::dh_da_time(const double * const s, 
                                                  const double * const alpha, 
                                                  double T,
                                                  double * const dhv) const
{
  for (size_t i = 0; i < nmodels(); i++) {
    size_t nhi = rules_[i]->nh();
    double * dhvi = new double [nhi*nhi];
    rules_[i]->dh_da_time(s, model_history_(alpha,i), T, dhvi);
    for (size_t a = 0; a < nhi; a++)
      for (size_t b = 0; b < nhi; b++)
        dhv[CINDEX(a+offsets_[i],b+offsets_[i],nh())] = dhvi[CINDEX(a,b,nhi)];

    delete [] dhvi;
  }
}

void SuperimposedViscoPlasticFlowRule::h_temp(const double * const s, 
                                              const double * const alpha,
                                              double T, double * const hv) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->h_temp(s, model_history_(alpha,i), T, model_history_(hv,i));
}

void SuperimposedViscoPlasticFlowRule::dh_ds_temp(const double * const s, 
                                                  const double * const alpha, 
                                                  double T, 
                                                  double * const dhv) const
{
  for (size_t i = 0; i < nmodels(); i++)
    rules_[i]->dh_ds_temp(s, model_history_(alpha,i), T, &(dhv[offsets_[i]*6]));
}

void SuperimposedViscoPlasticFlowRule::dh_da_temp(const double * const s, 
                                                  const double * const alpha,
                                                  double T,
                                                  double * const dhv) const
{
  for (size_t i = 0; i < nmodels(); i++) {
    size_t nhi = rules_[i]->nh();
    double * dhvi = new double [nhi*nhi];
    rules_[i]->dh_da_temp(s, model_history_(alpha,i), T, dhvi);
    for (size_t a = 0; a < nhi; a++)
      for (size_t b = 0; b < nhi; b++)
        dhv[CINDEX(a+offsets_[i],b+offsets_[i],nh())] = dhvi[CINDEX(a,b,nhi)];

    delete [] dhvi;
  }
}

double * SuperimposedViscoPlasticFlowRule::model_history_(double * const h, 
                                                          size_t i) const
{
  return &(h[offsets_[i]]);
}

const double * SuperimposedViscoPlasticFlowRule::model_history_(const double * const h, 
                                                                size_t i) const
{
  return &(h[offsets_[i]]);
}

GFlow::GFlow(ParameterSet & params) :
    NEMLObject(params)
{

}

// Various g(s) implementations
GPowerLaw::GPowerLaw(ParameterSet & params) :
    GFlow(params),
    n_(params.get_object_parameter<Interpolate>("n")), 
    eta_(params.get_object_parameter<Interpolate>("eta"))
{

}

std::string GPowerLaw::type()
{
  return "GPowerLaw";
}

ParameterSet GPowerLaw::parameters()
{
  ParameterSet pset(GPowerLaw::type());

  pset.add_parameter<NEMLObject>("n");
  pset.add_parameter<NEMLObject>("eta");

  return pset;
}

std::unique_ptr<NEMLObject> GPowerLaw::initialize(ParameterSet & params)
{
  return neml::make_unique<GPowerLaw>(params); 
}

double GPowerLaw::g(double f, double T) const
{
  return pow(f / eta_->value(T), n_->value(T));
}

double GPowerLaw::dg(double f, double T) const
{
  return n_->value(T) * pow(f / eta_->value(T), n_->value(T) - 1.0) / 
      eta_->value(T);
}

double GPowerLaw::n(double T) const
{
  return n_->value(T);
}

double GPowerLaw::eta(double T) const
{
  return eta_->value(T);
}

PerzynaFlowRule::PerzynaFlowRule(ParameterSet & params) :
    ViscoPlasticFlowRule(params),
    surface_(params.get_object_parameter<YieldSurface>("surface")), 
    hardening_(params.get_object_parameter<HardeningRule>("hardening")), 
    g_(params.get_object_parameter<GFlow>("g"))
{
  cache_history_(); 
}

std::string PerzynaFlowRule::type()
{
  return "PerzynaFlowRule";
}

ParameterSet PerzynaFlowRule::parameters()
{
  ParameterSet pset(PerzynaFlowRule::type());

  pset.add_parameter<NEMLObject>("surface");
  pset.add_parameter<NEMLObject>("hardening");
  pset.add_parameter<NEMLObject>("g");

  return pset;
}

std::unique_ptr<NEMLObject> PerzynaFlowRule::initialize(ParameterSet & params)
{
  return neml::make_unique<PerzynaFlowRule>(params); 
}

void PerzynaFlowRule::populate_hist(History & hist) const
{
  if (surface_->nhist() != hardening_->nhist()) {
    throw NEMLError("Hardening model and flow surface are not compatible");
  }
  
  // We want to make sure the variable prefix remains the same
  hardening_->set_variable_prefix(get_variable_prefix());

  hardening_->populate_hist(hist);
}

void PerzynaFlowRule::init_hist(History & hist) const
{
  hardening_->init_hist(hist);
}

// Rate rule
void PerzynaFlowRule::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  double fv;
  surface_->f(s, q, T, fv);

  if (fv > 0.0) {
    yv = g_->g(fabs(fv), T);
  }
  else {
    yv = 0.0;
  }

}

void PerzynaFlowRule::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  double fv;
  surface_->f(s, q, T, fv);

  std::fill(dyv, dyv + 6, 0.0);

  if (fv > 0.0) {
    double dgv = g_->dg(fabs(fv), T);
    surface_->df_ds(s, q, T, dyv);
    for (int i=0; i<6; i++) {
      dyv[i] *= dgv;
    }
  }
  
}

void PerzynaFlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  double fv;
  surface_->f(s, q, T, fv);

  std::fill(dyv, dyv + nh(), 0.0);

  if (fv > 0.0) {
    double dgv = g_->dg(fabs(fv), T);
    
    std::vector<double> jacv(nh()*nh());
    double * jac = &jacv[0];
    hardening_->dq_da(alpha, T, jac);
    
    std::vector<double> rdv(nh());
    double * rd = &rdv[0];
    surface_->df_dq(s, q, T, rd);

    mat_vec_trans(jac, nh(), rd, nh(), dyv);

    for (size_t i=0; i<nh(); i++) {
      dyv[i] *= dgv;
    }
  }


}

// Flow rule
void PerzynaFlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_ds(s, q, T, gv);
}

void PerzynaFlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_dsds(s, q, T, dgv);
}

void PerzynaFlowRule::dg_da(const double * const s, const double * const alpha, double T,
             double * const dgv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
  
  std::vector<double> jacv(nh() * nh());
  double * jac = &jacv[0];
  hardening_->dq_da(alpha, T, jac);
  
  std::vector<double> ddv(6*nh());
  double * dd = &ddv[0];
  surface_->df_dsdq(s, q, T, dd);

  mat_mat(6, nh(), nh(), dd, jac, dgv);
}

// Hardening rule
void PerzynaFlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_dq(s, q, T, hv);
}

void PerzynaFlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_dqds(s, q, T, dhv);
}

void PerzynaFlowRule::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::vector<double> qv(nh());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
  
  std::vector<double> jacv(nh() * nh());
  double * jac = &jacv[0];
  hardening_->dq_da(alpha, T, jac);
  
  std::vector<double> ddv(nh() * nh());
  double * dd = &ddv[0];
  surface_->df_dqdq(s, q, T, dd);

  mat_mat(nh(), nh(), nh(), dd, jac, dhv);
}

LinearViscousFlow::LinearViscousFlow(ParameterSet & params) :
    ViscoPlasticFlowRule(params),
    surface_(params.get_object_parameter<YieldSurface>("surface")),
    eta_(params.get_object_parameter<Interpolate>("eta"))
{
  cache_history_(); 
}

std::string LinearViscousFlow::type()
{
  return "LinearViscousFlow";
}

ParameterSet LinearViscousFlow::parameters()
{
  ParameterSet pset(LinearViscousFlow::type());

  pset.add_parameter<NEMLObject>("surface");
  pset.add_parameter<NEMLObject>("eta");

  return pset;
}

std::unique_ptr<NEMLObject> LinearViscousFlow::initialize(ParameterSet & params)
{
  return neml::make_unique<LinearViscousFlow>(params); 
}

void LinearViscousFlow::init_hist(History & h) const
{
  return;
}

void LinearViscousFlow::populate_hist(History & h) const
{
  return;
}

// Zero history of the correct size
std::vector<double> LinearViscousFlow::fake_hist_() const
{
  std::vector<double> h(surface_->nhist());
  std::fill(h.begin(), h.end(), 0.0);
  return h;
}

// Rate rule
void LinearViscousFlow::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  auto h = fake_hist_();

  double fv;
  surface_->f(s, &h[0], T, fv);
  
  yv = 3.0 / 2.0 * fv / eta_->value(T);
}

void LinearViscousFlow::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  auto h = fake_hist_();

  surface_->df_ds(s, &h[0], T, dyv);
  for (size_t i = 0; i < 6; i++)
    dyv[i] = 3.0 / 2.0 * dyv[i] / eta_->value(T);
}

void LinearViscousFlow::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  return;
}

// Flow rule
void LinearViscousFlow::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  auto h = fake_hist_();
  surface_->df_ds(s, &h[0], T, gv);
}

void LinearViscousFlow::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  auto h = fake_hist_();
  surface_->df_dsds(s, &h[0], T, dgv);
}

void LinearViscousFlow::dg_da(const double * const s, const double * const alpha, double T,
             double * const dgv) const
{
  return;
}

// Hardening rule
void LinearViscousFlow::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  return;
}

void LinearViscousFlow::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return;
}

void LinearViscousFlow::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  return;
}

FluidityModel::FluidityModel(ParameterSet & params) :
    NEMLObject(params)
{

}

// Begin Chaboche
ConstantFluidity::ConstantFluidity(ParameterSet & params) :
    FluidityModel(params),
    eta_(params.get_object_parameter<Interpolate>("eta"))
{

}

std::string ConstantFluidity::type()
{
  return "ConstantFluidity";
}

ParameterSet ConstantFluidity::parameters()
{
  ParameterSet pset(ConstantFluidity::type());

  pset.add_parameter<NEMLObject>("eta");

  return pset;
}

std::unique_ptr<NEMLObject> ConstantFluidity::initialize(ParameterSet & params)
{
  return neml::make_unique<ConstantFluidity>(params); 
}


double ConstantFluidity::eta(double a, double T) const
{
  return eta_->value(T);
}

double ConstantFluidity::deta(double a, double T) const
{
  return 0.0;
}

SaturatingFluidity::SaturatingFluidity(ParameterSet & params):
    FluidityModel(params),
    K0_(params.get_object_parameter<Interpolate>("K0")),
    A_(params.get_object_parameter<Interpolate>("A")),
    b_(params.get_object_parameter<Interpolate>("b"))
{

}

std::string SaturatingFluidity::type()
{
  return "SaturatingFluidity";
}

ParameterSet SaturatingFluidity::parameters()
{
  ParameterSet pset(SaturatingFluidity::type());

  pset.add_parameter<NEMLObject>("K0");
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("b");

  return pset;
}

std::unique_ptr<NEMLObject> SaturatingFluidity::initialize(ParameterSet & params)
{
  return neml::make_unique<SaturatingFluidity>(params); 
}


double SaturatingFluidity::eta(double a, double T) const
{
  double K0 = K0_->value(T);
  double A = A_->value(T);
  double b = b_->value(T);

  return K0 + A * (1.0 - exp(-b * a));
}

double SaturatingFluidity::deta(double a, double T) const
{
  double A = A_->value(T);
  double b = b_->value(T);

  return A * b * exp(-b * a);
}

ChabocheFlowRule::ChabocheFlowRule(ParameterSet & params) :
    ViscoPlasticFlowRule(params),
    surface_(params.get_object_parameter<YieldSurface>("surface")), 
    hardening_(params.get_object_parameter<NonAssociativeHardening>("hardening")), 
    fluidity_(params.get_object_parameter<FluidityModel>("fluidity")), 
    n_(params.get_object_parameter<Interpolate>("n")),
    prefactor_(params.get_object_parameter<Interpolate>("prefactor")), 
    recovery_(false)
{
  cache_history_();
}

std::string ChabocheFlowRule::type()
{
  return "ChabocheFlowRule";
}

ParameterSet ChabocheFlowRule::parameters() {
  ParameterSet pset(ChabocheFlowRule::type());

  pset.add_parameter<NEMLObject>("surface");
  pset.add_parameter<NEMLObject>("hardening");
  pset.add_parameter<NEMLObject>("fluidity");
  pset.add_parameter<NEMLObject>("n");
  pset.add_optional_parameter<NEMLObject>(
      "prefactor", make_constant(1.0));

  return pset;
}

std::unique_ptr<NEMLObject> ChabocheFlowRule::initialize(ParameterSet & params) {
  return neml::make_unique<ChabocheFlowRule>(params);
}

void ChabocheFlowRule::populate_hist(History & hist) const
{
  if (surface_->nhist() != hardening_->ninter()) {
    throw NEMLError("Hardening model and flow surface are not compatible");
  }
  
  // Make sure the variable prefix stays the same
  hardening_->set_variable_prefix(get_variable_prefix());

  hardening_->populate_hist(hist);
}

void ChabocheFlowRule::init_hist(History & hist) const
{
  hardening_->init_hist(hist);
}

// Rate rule
void ChabocheFlowRule::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
  
  double fv;
  surface_->f(s, q, T, fv);

  if (fv > 0.0) {
    double eta = sqrt(2.0 / 3.0) * fluidity_->eta(alpha[0], T);
    yv = sqrt(3.0 / 2.0) * pow(fv / eta, n_->value(T)) * prefactor_->value(T);
  } else {
    yv = 0.0;
  }

}

void ChabocheFlowRule::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
 
  double fv;
  surface_->f(s, q, T, fv);
  
  std::fill(dyv, dyv+6, 0.0);

  if (fv > 0.0) {
    surface_->df_ds(s, q, T, dyv);
    double eta = sqrt(2.0 / 3.0) * fluidity_->eta(alpha[0], T);
    double mv = sqrt(3.0 / 2.0) * pow(fv / eta, n_->value(T) - 1.0) *
                n_->value(T) / eta * prefactor_->value(T);
    for (int i = 0; i < 6; i++)
      dyv[i] *= mv;
  }

}

void ChabocheFlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  double fv;
  surface_->f(s, q, T, fv);

  std::fill(dyv, dyv + nh(), 0.0);

  if (fv > 0.0) {
    std::vector<double> jacv(hardening_->ninter() * nh());
    double * jac = &jacv[0];
    hardening_->dq_da(alpha, T, jac);
    
    std::vector<double> dqv(hardening_->ninter());
    double * dq = &dqv[0];
    surface_->df_dq(s, q, T, dq);

    mat_vec_trans(jac, nh(), dq, hardening_->ninter(), dyv);

    double eta = sqrt(2.0/3.0) * fluidity_->eta(alpha[0], T);
    double mv = sqrt(3.0 / 2.0) * pow(fv / eta, n_->value(T) - 1.0) *
                n_->value(T) / eta * prefactor_->value(T);
    for (size_t i = 0; i < nh(); i++)
      dyv[i] *= mv;

    double mv2 = -sqrt(3.0 / 2.0) * fv * pow(fv / eta, n_->value(T) - 1.0) *
                 n_->value(T) / (eta * eta) * prefactor_->value(T);

    double deta = sqrt(2.0/3.0) * fluidity_->deta(alpha[0], T);
    dyv[0] += deta * mv2;
  }


}

// Flow rule
void ChabocheFlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_ds(s, q, T, gv);
}

void ChabocheFlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);

  surface_->df_dsds(s, q, T, dgv);
}

void ChabocheFlowRule::dg_da(const double * const s, const double * const alpha, double T,
             double * const dgv) const
{
  std::vector<double> qv(hardening_->ninter());
  double * q = &qv[0];
  hardening_->q(alpha, T, q);
  
  std::vector<double> jacv(hardening_->ninter() * nh());
  double * jac = &jacv[0];
  hardening_->dq_da(alpha, T, jac);
  
  std::vector<double> ddv(6 * hardening_->ninter());
  double * dd = &ddv[0];
  surface_->df_dsdq(s, q, T, dd);

  mat_mat(6, nh(), hardening_->ninter(), dd, jac, dgv);
}

// Hardening rule
void ChabocheFlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  hardening_->h(s, alpha, T, hv);
}

void ChabocheFlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  hardening_->dh_ds(s, alpha, T, dhv);
}

void ChabocheFlowRule::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  hardening_->dh_da(s, alpha, T, dhv);
}

// Hardening rule wrt time
void ChabocheFlowRule::h_time(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  hardening_->h_time(s, alpha, T, hv);
}

void ChabocheFlowRule::dh_ds_time(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  hardening_->dh_ds_time(s, alpha, T, dhv);
}

void ChabocheFlowRule::dh_da_time(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  hardening_->dh_da_time(s, alpha, T, dhv);
}

// Hardening rule wrt temperature
void ChabocheFlowRule::h_temp(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  hardening_->h_temp(s, alpha, T, hv);
}

void ChabocheFlowRule::dh_ds_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  hardening_->dh_ds_temp(s, alpha, T, dhv);
}

void ChabocheFlowRule::dh_da_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  hardening_->dh_da_temp(s, alpha, T, dhv);
}

YaguchiGr91FlowRule::YaguchiGr91FlowRule(ParameterSet & params) : 
    ViscoPlasticFlowRule(params)
{
  cache_history_();
}

std::string YaguchiGr91FlowRule::type()
{
  return "YaguchiGr91FlowRule";
}

ParameterSet YaguchiGr91FlowRule::parameters()
{
  ParameterSet pset(YaguchiGr91FlowRule::type());

  return pset;
}

std::unique_ptr<NEMLObject> YaguchiGr91FlowRule::initialize(ParameterSet & params)
{
  return neml::make_unique<YaguchiGr91FlowRule>(params); 
}


void YaguchiGr91FlowRule::populate_hist(History & hist) const
{
  // Order:
  // 0-5  X1
  // 6-11 X2
  // 12   Q
  // 13   sa
  hist.add<Symmetric>(prefix("X1"));
  hist.add<Symmetric>(prefix("X2"));
  hist.add<double>(prefix("Q"));
  hist.add<double>(prefix("sa"));
}

void YaguchiGr91FlowRule::init_hist(History & hist) const
{
  hist.get<Symmetric>(prefix("X1")) = Symmetric::zero();
  hist.get<Symmetric>(prefix("X2")) = Symmetric::zero();
  hist.get<double>(prefix("Q")) = 0.0;
  hist.get<double>(prefix("sa")) = 0.0;
}

// Rate rule
void YaguchiGr91FlowRule::y(const double* const s, const double* const alpha, double T,
              double & yv) const
{
  double nT = n(T);
  double DT = D(T);
  double sa = alpha[13];

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  yv = (J2_(dS) - sa) / DT;
  if (yv > 0.0) {
    yv = pow(fabs(yv), nT);
  }
  else {
    yv = 0.0;
  }
}

void YaguchiGr91FlowRule::dy_ds(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  std::fill(dyv, dyv+6, 0.0);

  double yi;
  y(s, alpha, T, yi);
  double nT = n(T);
  double DT = D(T);
  double sa = alpha[13];

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);
  
  if (yi > 0.0) {
    double j2 = J2_(dS);
    double sp = (j2 - sa) / DT;
    sp = pow(fabs(sp), nT - 1.0) * nT * copysign(1.0, sp) / DT;
    dev_vec_deriv_(dS, dyv);
    for (int i=0; i<6; i++) {
      dyv[i] *= 3.0/2.0 / j2 * sp;
    }
  }
  else {
    std::fill(dyv, dyv+6, 0.0);
  }

}

void YaguchiGr91FlowRule::dy_da(const double* const s, const double* const alpha, double T,
              double * const dyv) const
{
  std::fill(dyv, dyv+nh(), 0.0);

  // General
  double yi;
  y(s, alpha, T, yi);
  double nT = n(T);
  double DT = D(T);
  double sa = alpha[13];

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);
  double j2 = J2_(dS);
  double q = (j2 - sa) / DT;

  if (yi > 0.0) {
    // Xs
    double sp = (j2 - sa) / DT;
    sp = pow(fabs(sp), nT - 1.0) * nT * copysign(1.0, sp) / DT;

    dev_vec_deriv_(dS, &dyv[0]);
    for (int i=0; i<6; i++) {
      dyv[i] *= -3.0/2.0 / j2 * sp;
    }

    dev_vec_deriv_(dS, &dyv[6]);
    for (int i=6; i<12; i++) {
      dyv[i] *= -3.0/2.0 / j2 * sp;
    }

    double dS2[6];
    sub_vec(s, &alpha[6], 6, dS2);

    // q (0)
    dyv[12] = 0.0;

    // sa
    dyv[13] = -nT * pow(fabs(q), nT - 1.0) * copysign(1.0, q) / DT;
  }
  else {
    std::fill(dyv, dyv+nh(), 0.0);
  }
  


}

// Flow rule
void YaguchiGr91FlowRule::g(const double * const s, const double * const alpha, double T,
              double * const gv) const
{
  std::fill(gv, gv+6, 0.0);

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  double Jn = J2_(dS);

  dev_vec(dS);

  if (Jn > 0.0) {
    for (int i=0; i<6; i++) {
      gv[i] = 3.0/2.0 * dS[i] / Jn;
    }
  }

}

void YaguchiGr91FlowRule::dg_ds(const double * const s, const double * const alpha, double T,
              double * const dgv) const
{
  std::fill(dgv, dgv+36, 0.0);

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      if (i==j) {
        dgv[CINDEX(i,j,6)] = 2.0/3.0;
      }
      else {
        dgv[CINDEX(i,j,6)] = -1.0/3.0;
      }
    }
  }
  for (int i=3; i<6; i++) {
    dgv[CINDEX(i,i,6)] = 1.0;
  }

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS); 
  double j2 = J2_(dS);
  double mdS[6];
  dev_vec_deriv_(dS, mdS);
  dev_vec(dS);
  
  for (int i=0; i<6; i++) {
    dS[i] *= 3.0 / (2.0 * pow(j2,2.0));
  }

  outer_update_minus(dS, 6, mdS, 6, dgv);

  for (int i=0; i<36; i++) {
    dgv[i] *= 3.0/(2.0 * j2);
  }

}

void YaguchiGr91FlowRule::dg_da(const double * const s, const double * const alpha, double T,
             double * const dgv) const
{
  // Only the X terms have derivatives
  std::fill(dgv, dgv+(6*nh()), 0.0);

  int nc = nh();

  // Bizarrely this is the easiest way to do this
  double deriv[36];
  dg_ds(s, alpha, T, deriv);

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      dgv[CINDEX(i,(j+0),nc)] = -deriv[CINDEX(i,j,6)];
      dgv[CINDEX(i,(j+6),nc)] = -deriv[CINDEX(i,j,6)];
    }
  }

}

// Hardening rule
void YaguchiGr91FlowRule::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  std::fill(hv, hv+nh(), 0.0);

  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  double Jn = J2_(dS);
  
  double n[6];
  dev_vec(dS);
  for (int i=0; i<6; i++) {
    n[i] = 3.0/2.0 * dS[i] / Jn;
  }

  double C1i = C1(T);
  double a1i = a10(T) - alpha[12];

  double C2i = C2(T);
  double a2i = a2(T);

  // X1
  for (int i=0; i<6; i++) {
    hv[i+0] = C1i * (2.0/3.0 * a1i * n[i] - alpha[i+0]);
  }
  
  // X2
  for (int i=0; i<6; i++) {
    hv[i+6] = C2i * (2.0/3.0 * a2i * n[i] - alpha[i+6]);
  }

  // Q
  double qi = q(T);
  double di = d(T);
  hv[12] = di*(qi - alpha[12]);

  // sa
  double bri = br(T);
  double bhi = bh(T);
  double Ai = A(T);
  double Bi = B(T);
  double yi;
  y(s, alpha, T, yi);

  if (fabs(yi) > log_tol_) {
    double sas = Ai + Bi * log10(yi);
    
    if (sas < 0.0) {
      sas = 0.0;
    }
    double bi;
    if ((sas - alpha[13]) >= 0.0) {
      bi = bhi;
    }
    else {
      bi = bri;
    }
    hv[13] = bi * (sas - alpha[13]);
  }
  else {
    hv[13] = 0.0;
  }

}

void YaguchiGr91FlowRule::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  // Only the X terms have derivatives
  std::fill(dhv, dhv + nh()*6, 0.0);

  double C1i = C1(T);
  double a1i = a10(T) - alpha[12];

  double C2i = C2(T);
  double a2i = a2(T);

  // Again, this is the easiest way to do this
  double deriv[36];
  dg_ds(s, alpha, T, deriv);

  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((i+0),j,6)] = deriv[CINDEX(i,j,6)] * 2.0/3.0 * C1i * a1i;
      dhv[CINDEX((i+6),j,6)] = deriv[CINDEX(i,j,6)] * 2.0/3.0 * C2i * a2i;
    }
  }

  // The derivative of the rate wrt to the stress goes into the last row
  double bri = br(T);
  double bhi = bh(T);
  double Ai = A(T);
  double Bi = B(T);
  double yi;
  y(s, alpha, T, yi);
  if (fabs(yi) > log_tol_) {
    double sas = Ai + Bi * log10(yi);
    if (sas > 0.0) {
      double bi;
      if ((sas - alpha[13]) >= 0.0) {
        bi = bhi;
      }
      else {
        bi = bri;
      }
      dy_ds(s, alpha, T, &dhv[CINDEX(13,0,6)]);
      for (int i=0; i<6; i++) {
        dhv[CINDEX(13,i,6)] = bi * Bi / (yi * log(10.0)) * dhv[CINDEX(13,i,6)];
      }
    }
  }

}

void YaguchiGr91FlowRule::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  // Fair number of cross-terms are zero
  int nhi = nh();
  std::fill(dhv, dhv+(nhi*nhi), 0.0);

  // Generic X terms
  std::vector<double> derivv(6*nhi);
  double * deriv = &derivv[0];
  dg_da(s, alpha, T, deriv);
  double C1i = C1(T);
  double a1i = a10(T) - alpha[12];

  double C2i = C2(T);
  double a2i = a2(T);

  for (int i=0; i<6; i++) {
    for (int j=0; j<nhi; j++) {
      dhv[CINDEX((i+0),j,nhi)] = 2.0/3.0 * a1i * deriv[CINDEX(i,j,nhi)];
      dhv[CINDEX((i+6),j,nhi)] = 2.0/3.0 * a2i * deriv[CINDEX(i,j,nhi)];
    }
  }

  for (int i=0; i<6; i++) {
    dhv[CINDEX((i+0),(i+0),nhi)] -= 1.0;
    dhv[CINDEX((i+6),(i+6),nhi)] -= 1.0;
  }


  for (int i=0; i<6; i++) {
    for (int j=0; j<nhi; j++) {
      dhv[CINDEX((i+0),j,nhi)] *= C1i;
      dhv[CINDEX((i+6),j,nhi)] *= C2i;
    }
  }

  // X1 has an extra term for the time-varying a1
  double X[6];
  std::fill(X, X+6, 0.0);
  add_vec(&alpha[0], &alpha[6], 6, X);
  double dS[6];
  sub_vec(s, X, 6, dS);

  double Jn = J2_(dS);
  
  double n[6];
  dev_vec(dS);
  for (int i=0; i<6; i++) {
    n[i] = 3.0/2.0 * dS[i] / Jn;
  }

  for (int i=0; i<6; i++) {
    dhv[CINDEX(i,12,nhi)] -= C1i*2.0/3.0*n[i];
  }

  // Q is nice and easy
  double di = d(T);
  dhv[CINDEX(12,12,nhi)] = -di;

  // There are two components to sa: the derivative of the rate wrt the history
  // and the derivative of the actual history term itself
  double bri = br(T);
  double bhi = bh(T);
  double Ai = A(T);
  double Bi = B(T);
  double yi;
  y(s, alpha, T, yi);
  
  if (fabs(yi) > log_tol_) {
    double sas = Ai + Bi * log10(yi);
    double bi;
    if ((sas - alpha[13]) >= 0.0) {
      bi = bhi;
    }
    else {
      bi = bri;
    }
    if (sas > 0.0) {
      dy_da(s, alpha, T, &dhv[CINDEX(13,0,nhi)]);
      for (int i=0; i<nhi; i++) {
        dhv[CINDEX(13,i,nhi)] = bi * Bi / (yi * log(10.0)) * dhv[CINDEX(13,i,nhi)];
      }
    }
    dhv[CINDEX(13,13,nhi)] += -bi;
  }

}

// Hardening rule wrt to time
void YaguchiGr91FlowRule::h_time(const double * const s, 
                                const double * const alpha, double T,
                                double * const hv) const
{
  std::fill(hv, hv+nh(), 0.0);
  
  double mi = m(T);

  // X1
  double g1i = g1(T);
  double J1 = J2_(&alpha[0]);
  for (int i=0; i<6; i++) {
    hv[i+0] = -g1i * pow(J1, mi-1.0) * alpha[i+0];
  }
  
  // X2
  double g2i = g2(T);
  double J2 = J2_(&alpha[6]);
  for (int i=0; i<6; i++) {
    hv[i+6] = -g2i * pow(J2, mi-1.0) * alpha[i+6];
  }

}

void YaguchiGr91FlowRule::dh_ds_time(const double * const s, 
                                    const double * const alpha, double T,
                                    double * const dhv) const
{
  // This is actually still zero
  std::fill(dhv, dhv+(nh()*6), 0.0);
}

void YaguchiGr91FlowRule::dh_da_time(const double * const s, 
                                    const double * const alpha, double T,
                                    double * const dhv) const
{
  int nhi = nh();
  std::fill(dhv, dhv+(nhi*nhi), 0.0);

  // This is non-zero
  double mi = m(T);

  // X1
  double g1i = g1(T);
  double J1 = J2_(&alpha[0]);
  double X1[6];
  std::copy(&alpha[0], &alpha[6], X1);
  double X1d[6];
  dev_vec_deriv_(X1, X1d);
  dev_vec(X1);
  for (int i=0; i<6; i++) {
    X1[i] *= (mi-1.0) * pow(J1,mi-3.0) * 3.0 / 2.0;
  }

  double dX1[36];
  std::fill(dX1, dX1+36, 0.0);
  for (int i=0; i<6; i++) {
    dX1[CINDEX(i,i,6)] = pow(J1, mi-1.0);
  }
  outer_update(X1, 6, X1d, 6, dX1);
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((i+0),(j+0),nhi)] = -g1i * dX1[CINDEX(i,j,6)];
    }
  }
  
  // X2
  double g2i = g2(T);
  double J2 = J2_(&alpha[6]);
  double X2[6];
  std::copy(&alpha[6], &alpha[12], X2);
  double X2d[6];
  dev_vec_deriv_(X2, X2d);
  dev_vec(X2);
  for (int i=0; i<6; i++) {
    X2[i] *= (mi-1.0) * pow(J2,mi-3.0) * 3.0 / 2.0;
  }

  double dX2[36];
  std::fill(dX2, dX2+36, 0.0);
  for (int i=0; i<6; i++) {
    dX2[CINDEX(i,i,6)] = pow(J2, mi-1.0);
  }
  outer_update(X2, 6, X2d, 6, dX2);
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      dhv[CINDEX((i+6),(j+6),nhi)] = -g2i * dX2[CINDEX(i,j,6)];
    }
  }
}


// Properties...

double YaguchiGr91FlowRule::D(double T) const
{
  return pow(-1.119e-9 + 5.145e-11*T - 5.450e-14*T*T, -1.0/2.85);
}

double YaguchiGr91FlowRule::n(double T) const
{
  return 2.850;
}

double YaguchiGr91FlowRule::a10(double T) const
{
 return 2.082e3 - 8.110*T + 1.321e-2 * T*T - 7.278e-6 * T*T*T;
}

double YaguchiGr91FlowRule::C2(double T) const
{
  if (T < 798.0) {
    return 2.0e2;
  }
  else {
    return -2.992e3 + 4.000 * T;
  }
}

double YaguchiGr91FlowRule::a2(double T) const
{
  if (T < 773.0) {
    return 1.100e2;
  }
  else {
    return 3.373e3 - 7.622 * T + 4.400e-3 * T*T;
  }
}

double YaguchiGr91FlowRule::g1(double T) const
{
  if (T < 773.0) {
    return 5.064e-137 * exp(3.545e-1 * T);
  }
  else {
    return 5.336e-33 * exp(4.470e-2 * T);
  }
}

double YaguchiGr91FlowRule::g2(double T) const
{
  if (T < 773.0) {
    return 8.572e-108 * exp(2.771e-1 * T);
  }
  else {
    return 2.817e-11 - 7.538e-14 * T + 5.039e-17 * T*T;
  }
}

double YaguchiGr91FlowRule::m(double T) const
{
  if (T < 673.0) {
    return 1.200e1;
  }
  else if ((673.0 <= T) && (T < 773)) {
    return 5.766e2 * exp(-5.754e-3 * T);
  }
  else {
    return 6.750;
  }
}

double YaguchiGr91FlowRule::br(double T) const
{
  return 1.0e3;
}

double YaguchiGr91FlowRule::bh(double T) const
{
  return 1.065e3 - 1.404 * T + 4.417e-4 * T*T;
}

double YaguchiGr91FlowRule::A(double T) const
{
  if (T < 673.0) {
    return -1.841e2 + 2.075e-1 * T;
  }
  else if ((673.0 <= T) && (T < 823)) {
    return -4.799e2 + 1.262 * T - 9.133e-4 * T*T;
  }
  else {
    return -6.000e1;
  }
}

double YaguchiGr91FlowRule::B(double T) const
{
  if (T < 773.0) {
    return -6.340e1 + 6.850e-2 * T;
  }
  else {
    return -1.730e1;
  }
}

double YaguchiGr91FlowRule::d(double T) const
{
  return 3.208 - 1.010e-2 * T + 1.012e-5 * T*T;
}

double YaguchiGr91FlowRule::q(double T) const
{
  if (T < 673.0) {
    return 9.500e1;
  }
  else if ((673.0 <= T) && (T < 823.0)) {
    return -2.467e2 + 7.320e-1 * T - 3.333e-4 * T*T;
  }
  else {
    return 1.300e2;
  }
}

double YaguchiGr91FlowRule::C1(double T) const
{
  if (T < 673.0) {
    return 1.500e3;
  }
  else if ((673.0 <= T) && (T < 773.0)) {
    return -2.879e4 + 45.0 * T;
  }
  else {
    return 6.000e3;
  }
}

// Couple of helpers
double YaguchiGr91FlowRule::J2_(const double * const v) const
{
  double cv[6];
  std::copy(v, v+6, cv);
  dev_vec(cv);
  return sqrt(3.0/2.0 * dot_vec(cv, cv, 6));
}

void YaguchiGr91FlowRule::dev_vec_deriv_(const double * const a, 
                                        double * const b) const
{
  for (int i=0; i<3; i++) {
    b[i] = 2.0/3.0 * a[i];
    for (int j=0; j<3; j++) {
      if (i == j) continue;
      b[i] -= a[j]/3.0;
    }
  }
  for (int i=3; i<6; i++) {
    b[i] = a[i];
  }

}

} // namespace neml
