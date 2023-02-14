#include "damage.h"
#include "elasticity.h"

#include <cmath>
#include <stdexcept>
#include <limits>

namespace neml {

NEMLDamagedModel_sd::NEMLDamagedModel_sd(ParameterSet & params) :
    NEMLModel_sd(params), 
    base_(params.get_object_parameter<NEMLModel_sd>("base"))
{

}

void NEMLDamagedModel_sd::populate_state(History & hist) const
{
  populate_damage(hist);
  base_->set_variable_prefix(get_variable_prefix());
  base_->populate_state(hist);
}

void NEMLDamagedModel_sd::init_state(History & hist) const
{
  init_damage(hist);
  base_->init_state(hist);
}

void NEMLDamagedModel_sd::set_elastic_model(std::shared_ptr<LinearElasticModel>
                                           emodel)
{
  elastic_ = emodel;
  base_->set_elastic_model(emodel);
}

NEMLScalarDamagedModel_sd::NEMLScalarDamagedModel_sd(ParameterSet & params) :
      NEMLDamagedModel_sd(params), 
      dmodel_(params.get_object_parameter<ScalarDamage>("damage")),
      rtol_(params.get_parameter<double>("rtol")), 
      atol_(params.get_parameter<double>("atol")), 
      miter_(params.get_parameter<int>("miter")),
      verbose_(params.get_parameter<bool>("verbose")), 
      linesearch_(params.get_parameter<bool>("linesearch")),
      ekill_(params.get_parameter<bool>("ekill")), 
      dkill_(params.get_parameter<double>("dkill")), 
      sfact_(params.get_parameter<double>("sfact"))
{
  cache_history_();
}

std::string NEMLScalarDamagedModel_sd::type()
{
  return "NEMLScalarDamagedModel_sd";
}

ParameterSet NEMLScalarDamagedModel_sd::parameters()
{
  ParameterSet pset(NEMLScalarDamagedModel_sd::type());
  
  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("base");
  pset.add_parameter<NEMLObject>("damage");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          make_constant(0.0));
  pset.add_optional_parameter<double>("rtol", 1.0e-6);
  pset.add_optional_parameter<double>("atol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<bool>("linesearch", false);
  pset.add_optional_parameter<bool>("truesdell", true);

  pset.add_optional_parameter<bool>("ekill", false);
  pset.add_optional_parameter<double>("dkill", 0.5);
  pset.add_optional_parameter<double>("sfact", 100000.0);

  return pset;
}

std::unique_ptr<NEMLObject> NEMLScalarDamagedModel_sd::initialize(ParameterSet & params)
{
  return neml::make_unique<NEMLScalarDamagedModel_sd>(params); 
}

void NEMLScalarDamagedModel_sd::update_sd_actual(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double & u_np1, double u_n,
    double & p_np1, double p_n)
{
  if (ekill_ and (h_n[0] >= dkill_)) {
    ekill_update_(T_np1, e_np1, s_np1, h_np1, h_n, A_np1, u_np1, u_n, p_np1, p_n);
    return;
  }

  // Make trial state
  SDTrialState tss;
  make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n, u_n, p_n, tss);
  
  // Call solve
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  solve(this, x, &tss, {rtol_, atol_, miter_, verbose_, linesearch_});
  
  // Do actual stress update
  double s_prime_np1[6];
  double s_prime_n[6];
  double A_prime_np1[36];
  std::copy(s_n, s_n+6, s_prime_n);
  for (int i=0; i<6; i++) s_prime_n[i] /= (1-h_n[0]);

  base_->update_sd_actual(e_np1, e_n, T_np1, T_n,
                   t_np1, t_n, s_prime_np1, s_prime_n,
                   &h_np1[1], &h_n[1],
                   A_prime_np1, u_np1, u_n, p_np1, p_n);


  for (int i=0; i<6; i++) s_np1[i] = (1-x[6]) * s_prime_np1[i];
  h_np1[0] = x[6];

  if (ekill_ and (h_np1[0] >= dkill_)) {
    ekill_update_(T_np1, e_np1, s_np1, h_np1, h_n, A_np1, u_np1, u_n, p_np1, p_n);
    return;
  }
  
  // Create the tangent
  tangent_(e_np1, e_n, s_np1, s_n,
                 T_np1, T_n, t_np1, t_n, 
                 x[6], h_n[0], A_prime_np1, A_np1);
}

size_t NEMLScalarDamagedModel_sd::ndamage() const
{
  return 1;
}

void NEMLScalarDamagedModel_sd::populate_damage(History & hist) const
{
  hist.add<double>(prefix("damage"));
}

void NEMLScalarDamagedModel_sd::init_damage(History & hist) const
{
  hist.get<double>(prefix("damage")) = dmodel_->d_init();
}

size_t NEMLScalarDamagedModel_sd::nparams() const
{
  return 7;
}

void NEMLScalarDamagedModel_sd::init_x(double * const x, TrialState * ts)
{
  SDTrialState * tss = static_cast<SDTrialState *>(ts);
  std::copy(tss->s_n, tss->s_n+6, x);
  // This provides a way for a particular model give a better guess at 
  // the initial damage for the first step.  Some models can be singular
  // for w = 0
  if (tss->w_n == 0.0) {
    x[6] = dmodel_->d_init();
  }
  else {
    x[6] = tss->w_n;
  }
}

void NEMLScalarDamagedModel_sd::RJ(const double * const x, TrialState * ts, 
                                  double * const R, double * const J)
{
  SDTrialState * tss = static_cast<SDTrialState *>(ts);
  const double * s_curr = x;
  double w_curr = x[6];
  double s_prime_curr[6];
  for (int i=0; i<6; i++)  s_prime_curr[i] = s_curr[i] / (1-w_curr);

  double s_prime_np1[6];
  double s_prime_n[6];
  double A_prime_np1[36];
  std::vector<double> h_np1_v(base_->nstate());
  double * h_np1 = &h_np1_v[0];
  double u_np1;
  double p_np1;
  
  std::copy(tss->s_n, tss->s_n+6, s_prime_n);
  for (int i=0; i<6; i++) s_prime_n[i] /= (1-tss->w_n);

  base_->update_sd_actual(tss->e_np1, tss->e_n, tss->T_np1, tss->T_n,
                   tss->t_np1, tss->t_n, s_prime_np1, s_prime_n,
                   h_np1, &tss->h_n[0],
                   A_prime_np1, u_np1, tss->u_n, p_np1, tss->p_n);
  
  for (int i=0; i<6; i++) R[i] = s_curr[i] - (1-w_curr) * s_prime_np1[i];

  double w_np1;
  dmodel_->damage(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n, tss->t_np1, tss->t_n, &w_np1);
  R[6] = w_curr - w_np1;

  std::fill(J, J+49, 0.0);
  for (int i=0; i<6; i++) {
    J[CINDEX(i,i,7)] = 1.0; 
  }
  for (int i=0; i<6; i++) {
    J[CINDEX(i,6,7)] = s_prime_np1[i];
  }

  double ws[6];
  dmodel_->ddamage_ds(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, ws);
  for (int i=0; i<6; i++) {
    J[CINDEX(6,i,7)] = -ws[i] / (1 - w_curr); 
  }
  
  double ww;
  dmodel_->ddamage_dd(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, &ww);
  
  J[CINDEX(6,6,7)] = 1.0 - ww - dot_vec(ws, s_curr, 6) / pow(1-w_curr,2.0);
}

void NEMLScalarDamagedModel_sd::make_trial_state(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n, double t_np1, double t_n,
    const double * const s_n, const double * const h_n,
    double u_n, double p_n,
    SDTrialState & tss)
{
  std::copy(e_np1, e_np1+6, tss.e_np1);
  std::copy(e_n, e_n+6, tss.e_n);
  tss.T_np1 = T_np1;
  tss.T_n = T_n;
  tss.t_np1 = t_np1;
  tss.t_n = t_n;
  std::copy(s_n, s_n+6, tss.s_n);
  tss.h_n.resize(base_->nstate());
  std::copy(h_n+1, h_n+base_->nstate()+1, tss.h_n.begin());
  tss.u_n = u_n;
  tss.p_n = p_n;
  tss.w_n = h_n[0];
}

void NEMLScalarDamagedModel_sd::tangent_(
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n, double t_np1, double t_n,
    double w_np1, double w_n, const double * const A_prime,
    double * const A)
{
  double s_prime_np1[6];
  for (int i=0; i<6; i++) s_prime_np1[i] = s_np1[i] / (1.0 - w_np1);
  double s_prime_n[6];
  for (int i=0; i<6; i++) s_prime_n[i] = s_n[i] / (1.0 - w_n);

  double dw_ds[6];
  dmodel_->ddamage_ds(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, dw_ds);
  double dw_de[6];
  dmodel_->ddamage_de(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, dw_de);
  double dw_dw;
  dmodel_->ddamage_dd(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, &dw_dw);

  double k1 = 1.0 - 1.0 / (1.0 - w_np1) * dot_vec(dw_ds, s_prime_np1, 6) - dw_dw;
  double B[36];
  std::fill(B, B+36, 0);
  for (int i=0; i<6; i++) B[CINDEX(i,i,6)] = 1.0;
  for (int i=0; i<6; i++) dw_ds[i] /= (k1 * (1.0 - w_np1));
  outer_update(s_prime_np1, 6, dw_ds, 6, B);
  invert_mat(B, 6);

  double C[36];
  std::copy(A_prime, A_prime+36, C);
  for (int i=0; i<36; i++) C[i] *= (1 - w_np1);
  for (int i=0; i<6; i++) dw_de[i] /= k1;
  outer_update_minus(s_prime_np1, 6, dw_de, 6, C);
  mat_mat(6, 6, 6, B, C, A);
}

void NEMLScalarDamagedModel_sd::ekill_update_(double T_np1, 
                                             const double * const e_np1,
                                             double * const s_np1,
                                             double * const h_np1,
                                             const double * const h_n,
                                             double * const A_np1, 
                                             double & u_np1, double u_n,
                                             double & p_np1, double p_n)
{
  std::copy(h_n, h_n + nstate(), h_np1);
  h_np1[0] = 1.0;
  elastic_->C(T_np1, A_np1);
  for (int i=0; i<36; i++) A_np1[i] /= sfact_;
  mat_vec(A_np1, 6, e_np1, 6, s_np1);
  if (u_n > 0.0) {
    p_np1 = p_n + u_n;
    u_np1 = 0.0;
  }
  else {
    p_np1 = p_n;
    u_np1 = u_n;
  }
}

double NEMLScalarDamagedModel_sd::get_damage(const double *const h_np1) 
{
  History h = gather_history_(h_np1);
  return h.get<double>(prefix("damage"));
}

bool NEMLScalarDamagedModel_sd::should_del_element(const double *const h_np1) 
{
  return get_damage(h_np1) > dkill_;
}

bool NEMLScalarDamagedModel_sd::is_damage_model() const 
{
  return true; 
}

ScalarDamage::ScalarDamage(ParameterSet & params) :
    NEMLObject(params),
    elastic_(params.get_object_parameter<LinearElasticModel>("elastic"))
{

}

ScalarDamageRate::ScalarDamageRate(ParameterSet & params) :
    ScalarDamage(params)
{

}

void ScalarDamageRate::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  damage_rate(d_np1, e_np1, s_np1, T_np1, t_np1, dd);

  *dd = (*dd) * (t_np1 - t_n) + d_n;
}

void ScalarDamageRate::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  ddamage_rate_dd(d_np1, e_np1, s_np1, T_np1, t_np1, dd);

  *dd = (*dd) * (t_np1 - t_n);
}

void ScalarDamageRate::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  ddamage_rate_de(d_np1, e_np1, s_np1, T_np1, t_np1, dd);
  
  double dt = t_np1 - t_n;
  for (size_t i = 0; i < 6; i++)
    dd[i] *= dt;
}

void ScalarDamageRate::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  ddamage_rate_ds(d_np1, e_np1, s_np1, T_np1, t_np1, dd);

  double dt = t_np1 - t_n;
  for (size_t i = 0; i < 6; i++)
    dd[i] *= dt;
}

CombinedDamage::CombinedDamage(ParameterSet & params) :
      ScalarDamage(params),
      models_(params.get_object_parameter_vector<ScalarDamage>("models"))
{
}

std::string CombinedDamage::type()
{
  return "CombinedDamage";
}

ParameterSet CombinedDamage::parameters()
{
  ParameterSet pset(CombinedDamage::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<std::vector<NEMLObject>>("models");

  return pset;
}

std::unique_ptr<NEMLObject> CombinedDamage::initialize(ParameterSet & params)
{
  return neml::make_unique<CombinedDamage>(params); 
}

void CombinedDamage::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  *dd = d_n;
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    double di;
    (*it)->damage(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, &di);
    *dd += (di - d_n);
  }
}

void CombinedDamage::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  *dd = 0.0;
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    double di;
    (*it)->ddamage_dd(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, &di);
    *dd += di;
  }
}

void CombinedDamage::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    double di[6];
    (*it)->ddamage_de(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, di);
    for (int i=0; i<6; i++) dd[i] += di[i];
  }
}

void CombinedDamage::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    double di[6];
    (*it)->ddamage_ds(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, di);
    for (int i=0; i<6; i++) dd[i] += di[i];
  }
}

ClassicalCreepDamage::ClassicalCreepDamage(ParameterSet & params) :
      ScalarDamageRate(params),
      A_(params.get_object_parameter<Interpolate>("A")),
      xi_(params.get_object_parameter<Interpolate>("xi")),
      phi_(params.get_object_parameter<Interpolate>("phi"))
{


}

std::string ClassicalCreepDamage::type()
{
  return "ClassicalCreepDamage";
}

ParameterSet ClassicalCreepDamage::parameters()
{
  ParameterSet pset(ClassicalCreepDamage::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("xi");
  pset.add_parameter<NEMLObject>("phi");

  return pset;
}

std::unique_ptr<NEMLObject> ClassicalCreepDamage::initialize(ParameterSet & params)
{
  return neml::make_unique<ClassicalCreepDamage>(params); 
}

void ClassicalCreepDamage::damage_rate(
    double d, const double * const e,
    const double * const s, double T, double t, 
    double * const dd) const
{
  double xi = xi_->value(T);
  double A = A_->value(T);
  double phi = phi_->value(T);

  double se = this->se(s);

  *dd = pow(se / A, xi) * pow(1.0 - d, -phi);
}

void ClassicalCreepDamage::ddamage_rate_dd(
    double d, const double * const e,
    const double * const s, double T, double t, 
    double * const dd) const
{
  double xi = xi_->value(T);
  double A = A_->value(T);
  double phi = phi_->value(T);

  double se = this->se(s);

  *dd = pow(se / A, xi) * phi * pow(1.0 - d, -(phi + 1.0));
}

void ClassicalCreepDamage::ddamage_rate_de(
    double d, const double * const e,
    const double * const s, double T, double t, 
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);
}

void ClassicalCreepDamage::ddamage_rate_ds(
    double d, const double * const e,
    const double * const s, double T, double t, 
    double * const dd) const
{
  double xi = xi_->value(T);
  double A = A_->value(T);
  double phi = phi_->value(T);

  double se = this->se(s);

  if (se == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return;
  }

  std::copy(s, s+6, dd);
  double s_m = (s[0] + s[1] + s[2]) / 3.0;
  for (int i=0; i<3; i++) dd[i] -= s_m;

  double sf = 3.0 * xi / (2.0 * A * se) * pow(se/A, xi - 1.0) * 
      pow(1 - d, -phi);
  for (int i=0; i<6; i++) dd[i] *= sf;
}

double ClassicalCreepDamage::se(const double * const s) const
{
  return sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) + 
               pow(s[2] - s[0], 2.0) + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + 
                                              pow(s[5], 2.0))) / 2.0);
}

EffectiveStress::EffectiveStress(ParameterSet & params) :
    NEMLObject(params)
{

}

VonMisesEffectiveStress::VonMisesEffectiveStress(ParameterSet & params) : 
    EffectiveStress(params)
{

}

std::string VonMisesEffectiveStress::type()
{
  return "VonMisesEffectiveStress";
}

ParameterSet VonMisesEffectiveStress::parameters()
{
  ParameterSet pset(VonMisesEffectiveStress::type());

  return pset;
}

std::unique_ptr<NEMLObject> VonMisesEffectiveStress::initialize(ParameterSet & params)
{
  return neml::make_unique<VonMisesEffectiveStress>(params);
}

void VonMisesEffectiveStress::effective(const double * const s, double & eff) const
{
  eff = sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) 
              + pow(s[2] - s[0], 2.0) 
              + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + pow(s[5], 2.0))) / 2.0);
}

void VonMisesEffectiveStress::deffective(const double * const s, double * const deff) const
{
  std::copy(s,s+6,deff);
  double mean = (s[0] + s[1] + s[2]) / 3.0;
  for (int i=0; i<3; i++) deff[i] -= mean;
  double se;
  effective(s, se);
  
  if (se == 0) return;

  for (int i=0; i<6; i++) deff[i] *= 3.0/2.0 / se;
}

MeanEffectiveStress::MeanEffectiveStress(ParameterSet & params) :
  EffectiveStress(params)
{
}

std::string MeanEffectiveStress::type() { return "MeanEffectiveStress"; }

ParameterSet MeanEffectiveStress::parameters() {
  ParameterSet pset(MeanEffectiveStress::type());

  return pset;
}

std::unique_ptr<NEMLObject>
MeanEffectiveStress::initialize(ParameterSet &params) {
  return neml::make_unique<MeanEffectiveStress>(params);
}

void MeanEffectiveStress::effective(const double *const s, double &eff) const {
  eff = (s[0] + s[1] + s[2]) / 3.0;
}

void MeanEffectiveStress::deffective(const double *const s,
                                    double *const deff) const {

  for (int i = 0; i < 6; i++) {
    if (i < 3)
      deff[i] = 1. / 3.;
    else
      deff[i] = 0;
  }
}


HuddlestonEffectiveStress::HuddlestonEffectiveStress(ParameterSet & params) :
    EffectiveStress(params),
    b_(params.get_parameter<double>("b"))
{

}

std::string HuddlestonEffectiveStress::type()
{
  return "HuddlestonEffectiveStress";
}

ParameterSet HuddlestonEffectiveStress::parameters()
{
  ParameterSet pset(HuddlestonEffectiveStress::type());

  pset.add_parameter<double>("b");

  return pset;
}

std::unique_ptr<NEMLObject> HuddlestonEffectiveStress::initialize(ParameterSet & params)
{
  return neml::make_unique<HuddlestonEffectiveStress>(params);
}

void HuddlestonEffectiveStress::effective(const double * const s, double & eff) const
{
  double sd[6];
  std::copy(s, s+6, sd);
  dev_vec(sd);

  double I1v = I1(s);
  double I2v = I2(s);
  double I2pv = I2(sd);
  double se = sqrt(-3.0 * I2pv);
  double ss = sqrt(-3.0 * I2pv + I2v);

  // check if ss = 0.0
  if ( ss == 0.)
  {
    eff = 0.;
    return;
  }

  eff = se * exp(b_*(I1v / ss - 1.0));
}

void HuddlestonEffectiveStress::deffective(const double * const s, double * const deff) const
{
  // Useful common factors
  double sd[6];
  std::copy(s, s+6, sd);
  dev_vec(sd);
  std::fill(deff, deff+6, 0.0);
  
  if (norm2_vec(s, 6) == 0.0) return;

  double I1v = I1(s);
  double I2v = I2(s);
  double I2pv = I2(sd);
  double se = sqrt(-3.0 * I2pv);
  double ss = sqrt(-3.0 * I2pv + I2v);

  // check if ss = 0.0
  if (ss == 0.0) return;
 
  double eff = se * exp(b_*(I1v / ss - 1.0));
  
  // I1 term
  double t1 = b_ * eff / sqrt(I2v - 3.0 * I2pv);
  for (int i=0; i<3; i++) {
    deff[i] += t1;
  }

  // I2 term
  double t2 = -b_ * eff * I1v / (2.0 * pow(I2v - 3.0 * I2pv, 3.0/2.0));
  for (int i=0; i<3; i++) {
    deff[i] += t2 * I1v;
  }
  for (int i=0; i<6; i++) {
    deff[i] -= t2 * s[i];
  }

  // I2p term
  double t3 = 0.5 * eff * (3.0*b_*I1v/pow(I2v - 3.0*I2pv, 3.0/2.0) + 1.0 / I2pv);
  for (int i=0; i<6; i++) {
    deff[i] -= t3 * sd[i];
  }
}

MaxPrincipalEffectiveStress::MaxPrincipalEffectiveStress(ParameterSet & params)
  :
      EffectiveStress(params)
{

}

std::string MaxPrincipalEffectiveStress::type()
{
  return "MaxPrincipalEffectiveStress";
}

ParameterSet MaxPrincipalEffectiveStress::parameters()
{
  ParameterSet pset(MaxPrincipalEffectiveStress::type());

  return pset;
}

std::unique_ptr<NEMLObject> MaxPrincipalEffectiveStress::initialize(ParameterSet & params)
{
  return neml::make_unique<MaxPrincipalEffectiveStress>(params);
}

void MaxPrincipalEffectiveStress::effective(const double * const s, double & eff) const
{
  double vals[3];

  eigenvalues_sym(s, vals);
  eff = vals[2];

  if (eff < 0.0) eff = 0.0;
}

void MaxPrincipalEffectiveStress::deffective(const double * const s, double * const deff) const
{
  double vectors[9];
  double values[3];
  eigenvalues_sym(s, values);

  if (values[2] < 0.0) {
    std::fill(deff, deff+6, 0.0);
    return;
  }

  eigenvectors_sym(s, vectors);
  double * v = &(vectors[6]);

  double full[9];
  double nf = 0.0;
  for (int i=0; i<3; i++) {
    nf += v[i] * v[i];
    for (int j=0; j<3; j++) {
      full[CINDEX(i,j,3)] = v[i]*v[j];
    }
  }

  if (nf != 0) {
    for (int i=0; i<9; i++) {
      full[i] /= nf;
    }
  }

  sym(full, deff);
}


MaxSeveralEffectiveStress::MaxSeveralEffectiveStress(ParameterSet & params) :
    EffectiveStress(params),
    measures_(params.get_object_parameter_vector<EffectiveStress>("measures"))
{

}

std::string MaxSeveralEffectiveStress::type()
{
  return "MaxSeveralEffectiveStress";
}

ParameterSet MaxSeveralEffectiveStress::parameters()
{
  ParameterSet pset(MaxSeveralEffectiveStress::type());

  pset.add_parameter<std::vector<NEMLObject>>("measures");

  return pset;
}

std::unique_ptr<NEMLObject> MaxSeveralEffectiveStress::initialize(ParameterSet & params)
{
  return neml::make_unique<MaxSeveralEffectiveStress>(params);
}

void MaxSeveralEffectiveStress::effective(const double * const s, double & eff) const
{
  size_t ind;
  select_(s, ind, eff);
}

void MaxSeveralEffectiveStress::deffective(const double * const s, double * const deff) const
{
  size_t ind;
  double eff;
  select_(s, ind, eff);

  measures_[ind]->deffective(s, deff);
}

void MaxSeveralEffectiveStress::select_(const double * const s, size_t & ind, double & value) const
{
  value = -std::numeric_limits<double>::infinity();
  ind = -1;

  double vi;

  for (size_t i = 0; i<measures_.size(); i++) {
    measures_[i]->effective(s, vi);
    if (vi > value) {
      value = vi;
      ind = i;
    }
  }
}

SumSeveralEffectiveStress::SumSeveralEffectiveStress(ParameterSet & params) :
    EffectiveStress(params),
    measures_(params.get_object_parameter_vector<EffectiveStress>("measures")),
    weights_(params.get_parameter<std::vector<double>>("weights"))
{
  if (measures_.size() != weights_.size()) {
    throw std::invalid_argument("Length of measures and weights must be the same");
  }
}

std::string SumSeveralEffectiveStress::type()
{
  return "SumSeveralEffectiveStress";
}

ParameterSet SumSeveralEffectiveStress::parameters()
{
  ParameterSet pset(SumSeveralEffectiveStress::type());

  pset.add_parameter<std::vector<NEMLObject>>("measures");
  pset.add_parameter<std::vector<double>>("weights");

  return pset;
}

std::unique_ptr<NEMLObject> SumSeveralEffectiveStress::initialize(ParameterSet & params)
{
  return neml::make_unique<SumSeveralEffectiveStress>(params);
}

void SumSeveralEffectiveStress::effective(const double * const s, double & eff) const
{
  double val;
  eff = 0.0;
  for (size_t i = 0; i < measures_.size(); i++) {
    measures_[i]->effective(s, val);
    eff += weights_[i] * val;
  }
}

void SumSeveralEffectiveStress::deffective(const double * const s, double * const deff) const
{
  double val[6];
  std::fill(deff, deff+6, 0.0);
  for (size_t i = 0; i < measures_.size(); i++) {
    measures_[i]->deffective(s, val);
    for (size_t j = 0; j < 6; j++) {
      deff[j] += weights_[i] * val[j];
    }
  }
}

ModularCreepDamage::ModularCreepDamage(ParameterSet & params) :
      ScalarDamageRate(params),
      A_(params.get_object_parameter<Interpolate>("A")), 
      xi_(params.get_object_parameter<Interpolate>("xi")),
      phi_(params.get_object_parameter<Interpolate>("phi")), 
      estress_(params.get_object_parameter<EffectiveStress>("estress"))
{

}

std::string ModularCreepDamage::type()
{
  return "ModularCreepDamage";
}

ParameterSet ModularCreepDamage::parameters()
{
  ParameterSet pset(ModularCreepDamage::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("xi");
  pset.add_parameter<NEMLObject>("phi");
  pset.add_parameter<NEMLObject>("estress");

  return pset;
}

std::unique_ptr<NEMLObject> ModularCreepDamage::initialize(ParameterSet & params)
{
  return neml::make_unique<ModularCreepDamage>(params); 
}

void ModularCreepDamage::damage_rate(
    double d, const double * const e,
    const double * const s,
    double T, double t,
    double * const dd) const
{
  double xi = xi_->value(T);
  double A = A_->value(T);
  double phi = phi_->value(T);

  double se;
  estress_->effective(s, se);
  
  *dd = pow(se / A, xi) * pow(1.0 - d, xi-phi);
}

void ModularCreepDamage::ddamage_rate_dd(
    double d, const double * const e,
    const double * const s,
    double T, double t,
    double * const dd) const
{
  double xi = xi_->value(T);
  double A = A_->value(T);
  double phi = phi_->value(T);

  double se;
  estress_->effective(s, se);

  *dd = pow(se / A, xi) * (phi-xi) * pow(1.0 - d, xi-phi - 1.0);
}

void ModularCreepDamage::ddamage_rate_de(
    double d, const double * const e,
    const double * const s,
    double T, double t,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);
}

void ModularCreepDamage::ddamage_rate_ds(
    double d, const double * const e,
    const double * const s,
    double T, double t,
    double * const dd) const
{
  double xi = xi_->value(T);
  double A = A_->value(T);
  double phi = phi_->value(T);

  double se;
  estress_->effective(s, se);

  if (se == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return;
  }

  double scalar = pow(se / A, xi - 1.0) * xi / A * pow(1.0 - d, xi-phi);
  
  estress_->deffective(s, dd);
  for (int i=0; i<6; i++) dd[i] *= scalar;
}

LarsonMillerCreepDamage::LarsonMillerCreepDamage(ParameterSet & params) :
      ScalarDamageRate(params),
      lmr_(params.get_object_parameter<LarsonMillerRelation>("lmr")), 
      estress_(params.get_object_parameter<EffectiveStress>("estress"))
{
}

std::string LarsonMillerCreepDamage::type()
{
  return "LarsonMillerCreepDamage";
}

ParameterSet LarsonMillerCreepDamage::parameters()
{
  ParameterSet pset(LarsonMillerCreepDamage::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("lmr");
  pset.add_parameter<NEMLObject>("estress");

  return pset;
}

std::unique_ptr<NEMLObject> LarsonMillerCreepDamage::initialize(ParameterSet & params)
{
  return neml::make_unique<LarsonMillerCreepDamage>(params); 
}

void LarsonMillerCreepDamage::damage_rate(
    double d, const double * const e,
    const double * const s,
    double T, double t,
    double * const dd) const
{
  double se;
  estress_->effective(s, se);

  // The usual problem
  if (se == 0.0) {
    *dd = 0.0;
    return;
  }

  double tR;
  lmr_->tR(se * (1-d), T, tR);
  
  *dd = 1.0 / tR;
}

void LarsonMillerCreepDamage::ddamage_rate_dd(
    double d, const double * const e,
    const double * const s,
    double T, double t,
    double * const dd) const
{
  double se;
  estress_->effective(s, se);

  // The usual problem
  if (se == 0.0) {
    *dd = 0.0;
    return;
  }

  double tR;
  lmr_->tR(se * (1-d), T, tR);

  double dtR;
  lmr_->dtR_ds(se * (1-d), T, dtR);

  *dd = (dtR * se) / (tR * tR);
}

void LarsonMillerCreepDamage::ddamage_rate_de(
    double d, const double * const e,
    const double * const s,
    double T, double t,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);
}

void LarsonMillerCreepDamage::ddamage_rate_ds(
    double d, const double * const e,
    const double * const s,
    double T, double t,
    double * const dd) const
{
  double se;
  estress_->effective(s, se);

  // The usual problem
  if (se == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return;
  }

  double tR;
  lmr_->tR(se * (1-d), T, tR);

  double dtR;
  lmr_->dtR_ds(se * (1-d), T, dtR);

  double dse = -dtR * (1-d) / (tR * tR);

  estress_->deffective(s, dd);
  for (int i=0; i<6; i++) dd[i] *= dse;
}

StandardScalarDamage::StandardScalarDamage(ParameterSet & params) :
      ScalarDamage(params) 
{

}

void StandardScalarDamage::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  f(s_np1, d_np1, T_np1, fval);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);
  *dd = d_n + fval * deps;
}

void StandardScalarDamage::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double df;
  df_dd(s_np1, d_np1, T_np1, df);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  *dd = df * deps;
}

void StandardScalarDamage::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  f(s_np1, d_np1, T_np1, fval);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  if (deps == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return;
  }
  
  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] - s_n[i];
    de[i] = e_np1[i] - e_n[i];
  }

  double S[36];
  elastic_->S(T_np1, S);

  double dee[36];
  mat_vec(S, 6, ds, 6, dee);

  for (int i=0; i<6; i++) {
    dd[i] = (2.0 * fval) / (3.0 * deps) * (de[i] - dee[i]); 
  }
}

void StandardScalarDamage::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  f(s_np1, d_np1, T_np1, fval);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  if (deps == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return;
  }

  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] - s_n[i];
    de[i] = e_np1[i] - e_n[i];
  }

  double S[36];
  elastic_->S(T_np1, S);

  double dee[36];
  mat_vec(S, 6, ds, 6, dee);

  double v1[6];
  for (int i=0; i<6; i++) {
    v1[i] = (2.0 * fval) / (3.0 * deps) * (dee[i] - de[i]); 
  }

  mat_vec(S, 6, v1, 6, dd);

  double dds[6];
  df_ds(s_np1, d_np1, T_np1, dds);

  for (int i=0; i<6; i++) {
    dd[i] = dds[i] * deps + dd[i];
  }
}

double StandardScalarDamage::dep(
    const double * const s_np1, const double * const s_n,
    const double * const e_np1, const double * const e_n,
    double T_np1) const
{
  double S[36];
  elastic_->S(T_np1, S);
  
  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] - s_n[i];
    de[i] = e_np1[i] - e_n[i];
  }

  double dee[36];
  mat_vec(S, 6, ds, 6, dee);

  double val_sq =  2.0/3.0 * (dot_vec(de, de, 6) + dot_vec(dee, dee, 6) -
                         2.0 * dot_vec(de, dee, 6));
  
  if (val_sq < 0.0) {
    return 0.0;
  }
  else {
    return sqrt(val_sq);
  }
}

WorkDamage::WorkDamage(ParameterSet & params) :
      ScalarDamage(params), 
      Wcrit_(params.get_object_parameter<Interpolate>("Wcrit")),
      n_(params.get_parameter<double>("n")),
      eps_(params.get_parameter<double>("eps")),
      log_(params.get_parameter<bool>("log"))
{
}

std::string WorkDamage::type()
{
  return "WorkDamage";
}

ParameterSet WorkDamage::parameters()
{
  ParameterSet pset(WorkDamage::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("Wcrit");
  pset.add_parameter<double>("n");
  pset.add_optional_parameter<double>("eps", 1e-30);
  pset.add_optional_parameter<bool>("log", false);

  return pset;
}

std::unique_ptr<NEMLObject> WorkDamage::initialize(ParameterSet & params)
{
  return neml::make_unique<WorkDamage>(params); 
}

void WorkDamage::damage(double d_np1, double d_n,
                   const double * const e_np1, const double * const e_n,
                   const double * const s_np1, const double * const s_n,
                   double T_np1, double T_n,
                   double t_np1, double t_n,
                   double * const dd) const
{
  if (d_np1 == 0.0) {
    *dd = d_n;
    return;
  }

  double wrate = workrate(e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n,
                          std::fabs(d_np1), d_n);
  if (wrate == 0.0) {
    *dd = d_n;
    return;
  }

  double dt = t_np1 - t_n;

  double val = Wcrit(wrate);

  *dd = d_n + n_ * std::pow(std::fabs(d_np1), (n_-1.0)/n_) * 
      wrate * dt / val;
}

void WorkDamage::ddamage_dd(double d_np1, double d_n,
                   const double * const e_np1, const double * const e_n,
                   const double * const s_np1, const double * const s_n,
                   double T_np1, double T_n,
                   double t_np1, double t_n,
                   double * const dd) const
{
  if (d_np1 == 0.0) {
    *dd = d_n;
    return;
  }

  double wrate = workrate(e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n,
                          std::fabs(d_np1), d_n);
  if (wrate == 0.0) {
    *dd = d_n;
    return;
  }

  double val = Wcrit(wrate);
  double deriv = dWcrit(wrate);
  double dt = t_np1 - t_n;

  double S[36];
  elastic_->S(T_np1, S);
  double e[6];
  mat_vec(S, 6, s_np1, 6, e);
  double f = dot_vec(s_np1, e, 6);
  double x = -wrate / (1.0 - std::fabs(d_np1)) + 
      (1.0 - std::fabs(d_np1)) * f / dt;

  double other = n_*std::pow(fabs(d_np1), (n_-1.0)/n_) * 
      dt / val * (1.0 - wrate / val *deriv) * x;
  *dd = (n_-1.0)*std::pow(std::fabs(d_np1), -1.0/n_) *
      wrate * dt / val + other;
}

void WorkDamage::ddamage_de(double d_np1, double d_n,
                   const double * const e_np1, const double * const e_n,
                   const double * const s_np1, const double * const s_n,
                   double T_np1, double T_n,
                   double t_np1, double t_n,
                   double * const dd) const
{
  double wrate = workrate(e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n,
                          std::fabs(d_np1), d_n);

  // Provide a sensible answer if Newton's method gives us a garbage 
  // value of the history variable.
  if ((d_np1 <= 0.0) || (wrate == 0.0)) {
    std::fill(dd, dd+6, 0.0);
    return;
  }

  double val = Wcrit(wrate);
  double dval = dWcrit(wrate);

  double fact = n_ * std::pow(std::fabs(d_np1), (n_-1.0)/n_) / val * 
      (1.0 - wrate / val * dval) * (1.0 - std::fabs(d_np1));

  for (size_t i = 0; i < 6; i++) {
    dd[i] = fact * s_np1[i];
  }
}

void WorkDamage::ddamage_ds(double d_np1, double d_n,
                   const double * const e_np1, const double * const e_n,
                   const double * const s_np1, const double * const s_n,
                   double T_np1, double T_n,
                   double t_np1, double t_n,
                   double * const dd) const
{
  double wrate = workrate(e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n,
                          std::fabs(d_np1), d_n);

  if ((d_np1 <= 0.0) || (wrate == 0.0)) {
    std::fill(dd, dd+6, 0.0);
    return;
  }

  double val = Wcrit(wrate);
  double dval = dWcrit(wrate);

  double S[36];
  elastic_->S(T_np1, S);
  
  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] * (1.0-d_np1) - s_n[i] * (1.0-d_n);
    de[i] = e_np1[i] - e_n[i];
  }

  double dee[6];
  mat_vec(S, 6, ds, 6, dee);

  double e[6];
  mat_vec(S, 6, s_np1, 6, e);

  double fact = n_ * std::pow(fabs(d_np1), (n_-1.0)/n_) / val * 
      (1.0 - wrate / val * dval)*(1.0-std::fabs(d_np1));

  for (size_t i = 0; i < 6; i++) {
    dd[i] = fact * (de[i] - dee[i] - (1.0-std::fabs(d_np1))*e[i]);
  }
}

double WorkDamage::workrate(
    const double * const strain_np1, const double * const strain_n,
    const double * const stress_np1, const double * const stress_n,
    double T_np1, double T_n, double t_np1, double t_n,
    double d_np1, double d_n) const
{
  double dt = t_np1 - t_n;
  if (dt <= 0.0) {
    return 0.0;
  }

  double S[36];
  elastic_->S(T_np1, S);

  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = stress_np1[i] *(1.0-d_np1) - stress_n[i]*(1.0-d_n);
    de[i] = strain_np1[i] - strain_n[i];
  }

  double dee[6];
  mat_vec(S, 6, ds, 6, dee);

  double dp[6];
  for (int i = 0; i < 6; i++) {
    dp[i] = de[i] - dee[i];
  }

  return fabs(dot_vec(stress_np1, dp, 6) / dt * (1.0 - d_np1));
}

double WorkDamage::Wcrit(double Wdot) const
{
  if (log_)
    return std::pow(10.0, Wcrit_->value(std::log10(Wdot)));

  return Wcrit_->value(Wdot);
}

double WorkDamage::dWcrit(double Wdot) const
{
  if (log_)
    return std::pow(10.0, Wcrit_->value(std::log10(Wdot))) * 
        Wcrit_->derivative(std::log10(Wdot)) / Wdot;

  return Wcrit_->derivative(Wdot);
}

PowerLawDamage::PowerLawDamage(ParameterSet & params) :
      StandardScalarDamage(params), 
      A_(params.get_object_parameter<Interpolate>("A")),
      a_(params.get_object_parameter<Interpolate>("a"))
{

}

std::string PowerLawDamage::type()
{
  return "PowerLawDamage";
}

ParameterSet PowerLawDamage::parameters()
{
  ParameterSet pset(PowerLawDamage::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("a");

  return pset;
}

std::unique_ptr<NEMLObject> PowerLawDamage::initialize(ParameterSet & params)
{
  return neml::make_unique<PowerLawDamage>(params); 
}

void PowerLawDamage::f(const double * const s_np1, double d_np1,
                                double T_np1, double & f) const
{
  double sev = se(s_np1);
  double A = A_->value(T_np1);
  double a = a_->value(T_np1);

  f = A * pow(sev, a);
}

void PowerLawDamage::df_ds(const double * const s_np1, double d_np1, double T_np1,
                                 double * const df) const
{
  double sev = se(s_np1);
  double A = A_->value(T_np1);
  double a = a_->value(T_np1);

  if (sev == 0.0) {
    std::fill(df, df+6, 0.0);
    return;
  }

  std::copy(s_np1, s_np1+6, df);
  double sm = (s_np1[0] + s_np1[1] + s_np1[2]) / 3.0;
  for (int i=0; i<3; i++) df[i] -= sm;
  for (int i=0; i<6; i++) df[i] *= (3.0 * A * a / 2.0 * pow(sev, a - 2.0));
}

void PowerLawDamage::df_dd(const double * const s_np1, double d_np1, double T_np1,
                                 double & df) const
{
  df = 0.0;
}

double PowerLawDamage::se(const double * const s) const
{
  return sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) + 
               pow(s[2] - s[0], 2.0) + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + 
                                              pow(s[5], 2.0))) / 2.0);
}

ExponentialWorkDamage::ExponentialWorkDamage(ParameterSet
                                                                       & params) :
      StandardScalarDamage(params), 
      W0_(params.get_object_parameter<Interpolate>("W0")),
      k0_(params.get_object_parameter<Interpolate>("k0")),
      af_(params.get_object_parameter<Interpolate>("af"))
{

}

std::string ExponentialWorkDamage::type()
{
  return "ExponentialWorkDamage";
}

ParameterSet ExponentialWorkDamage::parameters()
{
  ParameterSet pset(ExponentialWorkDamage::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("W0");
  pset.add_parameter<NEMLObject>("k0");
  pset.add_parameter<NEMLObject>("af");

  return pset;
}

std::unique_ptr<NEMLObject> ExponentialWorkDamage::initialize(ParameterSet & params)
{
  return neml::make_unique<ExponentialWorkDamage>(params); 
}

void ExponentialWorkDamage::f(const double * const s_np1, double d_np1,
                                double T_np1, double & f) const
{
  double sev = se(s_np1);
  double W0 = W0_->value(T_np1);
  double k0 = k0_->value(T_np1);
  double af = af_->value(T_np1);
  
  // Sign, odd case during iteration
  if ((d_np1 + k0) < 0.0) {
    f = 0.0;
  }
  else {
    f = pow(d_np1 + k0, af) / W0 * sev;
  }
}

void ExponentialWorkDamage::df_ds(const double * const s_np1, double d_np1, double T_np1,
                                 double * const df) const
{
  double sev = se(s_np1);
  double W0 = W0_->value(T_np1);
  double k0 = k0_->value(T_np1);
  double af = af_->value(T_np1);

  if (sev == 0.0) {
    std::fill(df, df+6, 0.0);
    return;
  }

  if ((d_np1 + k0) < 0.0) {
    std::fill(df, df+6, 0.0);
    return;
  }

  std::copy(s_np1, s_np1+6, df);
  double sm = (s_np1[0] + s_np1[1] + s_np1[2]) / 3.0;
  for (int i=0; i<3; i++) df[i] -= sm;
  for (int i=0; i<6; i++) df[i] *= (3.0 * pow(d_np1 + k0, af) / (2.0 * sev * W0) );
}

void ExponentialWorkDamage::df_dd(const double * const s_np1, double d_np1, double T_np1,
                                 double & df) const
{
  double sev = se(s_np1);
  double W0 = W0_->value(T_np1);
  double k0 = k0_->value(T_np1);
  double af = af_->value(T_np1);

  if ((d_np1 + k0) < 0.0) {
    df = 0.0;
    return;
  }

  df = af * pow(d_np1 + k0, af - 1.0) * sev / W0;
}

double ExponentialWorkDamage::se(const double * const s) const
{
  return sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) + 
               pow(s[2] - s[0], 2.0) + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + 
                                              pow(s[5], 2.0))) / 2.0);
}

}
