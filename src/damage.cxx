#include "damage.h"
#include "elasticity.h"

#include <cmath>
#include <limits>

namespace neml {

NEMLDamagedModel_sd::NEMLDamagedModel_sd(ParameterSet & params) :
    NEMLModel_sd(params), 
    base_(params.get_object_parameter<NEMLModel_sd>("base"))
{

}

size_t NEMLDamagedModel_sd::nhist() const
{
  return ndamage() + base_->nhist();
}

int NEMLDamagedModel_sd::init_hist(double * const hist) const
{
  int res = init_damage(hist);
  if (res != 0) return res;
  return base_->init_hist(&hist[ndamage()]);
}

int NEMLDamagedModel_sd::set_elastic_model(std::shared_ptr<LinearElasticModel>
                                           emodel)
{
  elastic_ = emodel;
  return base_->set_elastic_model(emodel);
}

NEMLScalarDamagedModel_sd::NEMLScalarDamagedModel_sd(ParameterSet & params) :
      NEMLDamagedModel_sd(params), 
      rtol_(params.get_parameter<double>("rtol")), 
      atol_(params.get_parameter<double>("atol")), 
      miter_(params.get_parameter<int>("miter")),
      verbose_(params.get_parameter<bool>("verbose")), 
      linesearch_(params.get_parameter<bool>("linesearch")),
      ekill_(params.get_parameter<bool>("ekill")), 
      dkill_(params.get_parameter<double>("dkill")), 
      sfact_(params.get_parameter<double>("sfact"))
{
}

int NEMLScalarDamagedModel_sd::update_sd(
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
    return ekill_update_(T_np1, e_np1, s_np1, h_np1, h_n, A_np1, u_np1, u_n, p_np1, p_n);
  }

  // Make trial state
  SDTrialState tss;
  int ier = make_trial_state(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n, u_n, p_n, tss);
  if (ier != SUCCESS) return ier;
  
  // Call solve
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  ier = solve(this, x, &tss, {rtol_, atol_, miter_, verbose_, linesearch_});
  if (ier != SUCCESS) return ier;
  
  // Do actual stress update
  double s_prime_np1[6];
  double s_prime_n[6];
  double A_prime_np1[36];
  std::copy(s_n, s_n+6, s_prime_n);
  for (int i=0; i<6; i++) s_prime_n[i] /= (1-h_n[0]);

  ier = base_->update_sd(e_np1, e_n, T_np1, T_n,
                   t_np1, t_n, s_prime_np1, s_prime_n,
                   &h_np1[1], &h_n[1],
                   A_prime_np1, u_np1, u_n, p_np1, p_n);
  if (ier != SUCCESS) return ier;


  for (int i=0; i<6; i++) s_np1[i] = (1-x[6]) * s_prime_np1[i];
  h_np1[0] = x[6];

  if (ekill_ and (h_np1[0] >= dkill_)) {
    return ekill_update_(T_np1, e_np1, s_np1, h_np1, h_n, A_np1, u_np1, u_n, p_np1, p_n);
  }
  
  // Create the tangent
  ier = tangent_(e_np1, e_n, s_np1, s_n,
                 T_np1, T_n, t_np1, t_n, 
                 x[6], h_n[0], A_prime_np1, A_np1);
  if (ier != SUCCESS) return ier;

  return 0;
}

size_t NEMLScalarDamagedModel_sd::ndamage() const
{
  return 1;
}

int NEMLScalarDamagedModel_sd::init_damage(double * const damage) const
{
  damage[0] = 0.0;
  return 0;
}

size_t NEMLScalarDamagedModel_sd::nparams() const
{
  return 7;
}

int NEMLScalarDamagedModel_sd::init_x(double * const x, TrialState * ts)
{
  SDTrialState * tss = static_cast<SDTrialState *>(ts);
  std::copy(tss->s_n, tss->s_n+6, x);
  // This provides a way for a particular model give a better guess at 
  // the initial damage for the first step.  Some models can be singular
  // for w = 0
  if (tss->w_n == 0.0) {
    x[6] = d_guess();
  }
  else {
    x[6] = tss->w_n;
  }

  return 0;
}

int NEMLScalarDamagedModel_sd::RJ(const double * const x, TrialState * ts, 
                                  double * const R, double * const J)
{
  SDTrialState * tss = static_cast<SDTrialState *>(ts);
  const double * s_curr = x;
  double w_curr = x[6];
  double s_prime_curr[6];
  for (int i=0; i<6; i++)  s_prime_curr[i] = s_curr[i] / (1-w_curr);

  int res;
  double s_prime_np1[6];
  double s_prime_n[6];
  double A_prime_np1[36];
  std::vector<double> h_np1_v(base_->nhist());
  double * h_np1 = &h_np1_v[0];
  double u_np1;
  double p_np1;
  
  std::copy(tss->s_n, tss->s_n+6, s_prime_n);
  for (int i=0; i<6; i++) s_prime_n[i] /= (1-tss->w_n);

  res = base_->update_sd(tss->e_np1, tss->e_n, tss->T_np1, tss->T_n,
                   tss->t_np1, tss->t_n, s_prime_np1, s_prime_n,
                   h_np1, &tss->h_n[0],
                   A_prime_np1, u_np1, tss->u_n, p_np1, tss->p_n);
  if (res != SUCCESS) return res;
  
  for (int i=0; i<6; i++) R[i] = s_curr[i] - (1-w_curr) * s_prime_np1[i];

  double w_np1;
  res = damage(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n, tss->t_np1, tss->t_n, &w_np1);
  if (res != SUCCESS) return res;
  R[6] = w_curr - w_np1;

  std::fill(J, J+49, 0.0);
  for (int i=0; i<6; i++) {
    J[CINDEX(i,i,7)] = 1.0; 
  }
  for (int i=0; i<6; i++) {
    J[CINDEX(i,6,7)] = s_prime_np1[i];
  }

  double ws[6];
  res = ddamage_ds(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, ws);
  if (res != SUCCESS) return res;
  for (int i=0; i<6; i++) {
    J[CINDEX(6,i,7)] = -ws[i] / (1 - w_curr); 
  }
  
  double ww;
  res = ddamage_dd(w_curr, tss->w_n, tss->e_np1, tss->e_n, s_prime_curr, s_prime_n,
         tss->T_np1, tss->T_n,
         tss->t_np1, tss->t_n, &ww);
  if (res != SUCCESS) return res;
  
  J[CINDEX(6,6,7)] = 1.0 - ww - dot_vec(ws, s_curr, 6) / pow(1-w_curr,2.0);

  return 0;
}

int NEMLScalarDamagedModel_sd::make_trial_state(
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
  tss.h_n.resize(base_->nhist());
  std::copy(h_n+1, h_n+base_->nhist()+1, tss.h_n.begin());
  tss.u_n = u_n;
  tss.p_n = p_n;
  tss.w_n = h_n[0];

  return 0;
}

int NEMLScalarDamagedModel_sd::tangent_(
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
  int ier = ddamage_ds(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, dw_ds);
  if (ier != SUCCESS) return ier;
  double dw_de[6];
  ier = ddamage_de(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, dw_de);
  if (ier != SUCCESS) return ier;
  double dw_dw;
  ier = ddamage_dd(w_np1, w_n, e_np1, e_n, s_prime_np1, s_prime_n, T_np1, T_n,
             t_np1, t_n, &dw_dw);
  if (ier != SUCCESS) return ier;

  double k1 = 1.0 - 1.0 / (1.0 - w_np1) * dot_vec(dw_ds, s_prime_np1, 6) - dw_dw;
  double B[36];
  std::fill(B, B+36, 0);
  for (int i=0; i<6; i++) B[CINDEX(i,i,6)] = 1.0;
  for (int i=0; i<6; i++) dw_ds[i] /= (k1 * (1.0 - w_np1));
  outer_update(s_prime_np1, 6, dw_ds, 6, B);
  ier = invert_mat(B, 6);
  if (ier != SUCCESS) return ier;

  double C[36];
  std::copy(A_prime, A_prime+36, C);
  for (int i=0; i<36; i++) C[i] *= (1 - w_np1);
  for (int i=0; i<6; i++) dw_de[i] /= k1;
  outer_update_minus(s_prime_np1, 6, dw_de, 6, C);
  mat_mat(6, 6, 6, B, C, A);

  return 0;
}

int NEMLScalarDamagedModel_sd::ekill_update_(double T_np1, 
                                             const double * const e_np1,
                                             double * const s_np1,
                                             double * const h_np1,
                                             const double * const h_n,
                                             double * const A_np1, 
                                             double & u_np1, double u_n,
                                             double & p_np1, double p_n)
{
  std::copy(h_n, h_n + nhist(), h_np1);
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
  return 0;
}

CombinedDamageModel_sd::CombinedDamageModel_sd(ParameterSet & params) :
      NEMLScalarDamagedModel_sd(params),
      models_(params.get_object_parameter_vector<NEMLScalarDamagedModel_sd>("models"))
{
}

std::string CombinedDamageModel_sd::type()
{
  return "CombinedDamageModel_sd";
}

ParameterSet CombinedDamageModel_sd::parameters()
{
  ParameterSet pset(CombinedDamageModel_sd::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<std::vector<NEMLObject>>("models");
  pset.add_parameter<NEMLObject>("base");

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

std::unique_ptr<NEMLObject> CombinedDamageModel_sd::initialize(ParameterSet & params)
{
  return neml::make_unique<CombinedDamageModel_sd>(params); 
}

int CombinedDamageModel_sd::damage(
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
    int ier = (*it)->damage(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, &di);
    if (ier != SUCCESS) return ier;
    *dd += (di - d_n);
  }
  return 0;
}

int CombinedDamageModel_sd::ddamage_dd(
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
    int ier = (*it)->ddamage_dd(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, &di);
    if (ier != SUCCESS) return ier;
    *dd += di;
  }

  return 0;
}

int CombinedDamageModel_sd::ddamage_de(
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
    int ier = (*it)->ddamage_de(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, di);
    if (ier != SUCCESS) return ier;
    for (int i=0; i<6; i++) dd[i] += di[i];
  }

  return 0;
}

int CombinedDamageModel_sd::ddamage_ds(
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
    int ier = (*it)->ddamage_ds(d_np1, d_n, e_np1, e_n, s_np1, s_n, T_np1, T_n,
               t_np1, t_n, di);
    if (ier != SUCCESS) return ier;
    for (int i=0; i<6; i++) dd[i] += di[i];
  }

  return 0;
}

int CombinedDamageModel_sd::set_elastic_model(std::shared_ptr<LinearElasticModel>
                                          emodel)
{
  elastic_ = emodel;
  base_->set_elastic_model(emodel);
  for (auto it = models_.begin(); it != models_.end(); ++it) {
    int ier = (*it)->set_elastic_model(emodel);
    if (ier != SUCCESS) return ier;
  }
  return 0;
}

ClassicalCreepDamageModel_sd::ClassicalCreepDamageModel_sd(ParameterSet & params) :
      NEMLScalarDamagedModel_sd(params),
      A_(params.get_object_parameter<Interpolate>("A")),
      xi_(params.get_object_parameter<Interpolate>("xi")),
      phi_(params.get_object_parameter<Interpolate>("phi"))
{


}

std::string ClassicalCreepDamageModel_sd::type()
{
  return "ClassicalCreepDamageModel_sd";
}

ParameterSet ClassicalCreepDamageModel_sd::parameters()
{
  ParameterSet pset(ClassicalCreepDamageModel_sd::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("xi");
  pset.add_parameter<NEMLObject>("phi");
  pset.add_parameter<NEMLObject>("base");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          make_constant(0.0));
  pset.add_optional_parameter<double>("rtol", 1.0e-8);
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

std::unique_ptr<NEMLObject> ClassicalCreepDamageModel_sd::initialize(ParameterSet & params)
{
  return neml::make_unique<ClassicalCreepDamageModel_sd>(params); 
}

int ClassicalCreepDamageModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se = this->se(s_np1);
  double dt = t_np1 - t_n;

  *dd = d_n + pow(se / A, xi) * pow(1.0 - d_np1, -phi) * dt;

  return 0;
}

int ClassicalCreepDamageModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se = this->se(s_np1);
  double dt = t_np1 - t_n;

  *dd = pow(se / A, xi) * phi * pow(1.0 - d_np1, -(phi + 1.0)) * dt;

  return 0;
}

int ClassicalCreepDamageModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);

  return 0;
}

int ClassicalCreepDamageModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se = this->se(s_np1);
  double dt = t_np1 - t_n;

  if (se == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  std::copy(s_np1, s_np1+6, dd);
  double s_m = (s_np1[0] + s_np1[1] + s_np1[2]) / 3.0;
  for (int i=0; i<3; i++) dd[i] -= s_m;

  double sf = 3.0 * xi / (2.0 * A * se) * pow(se/A, xi - 1.0) * 
      pow(1 - d_np1, -phi) * dt;
  for (int i=0; i<6; i++) dd[i] *= sf;
  
  return 0;
}

double ClassicalCreepDamageModel_sd::se(const double * const s) const
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

int VonMisesEffectiveStress::effective(const double * const s, double & eff) const
{
  eff = sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) 
              + pow(s[2] - s[0], 2.0) 
              + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + pow(s[5], 2.0))) / 2.0);

  return 0;
}

int VonMisesEffectiveStress::deffective(const double * const s, double * const deff) const
{
  std::copy(s,s+6,deff);
  double mean = (s[0] + s[1] + s[2]) / 3.0;
  for (int i=0; i<3; i++) deff[i] -= mean;
  double se;
  effective(s, se);
  
  if (se == 0) return 0;

  for (int i=0; i<6; i++) deff[i] *= 3.0/2.0 / se;

  return 0;
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

int MeanEffectiveStress::effective(const double *const s, double &eff) const {
  eff = (s[0] + s[1] + s[2]) / 3.0;

  return 0;
}

int MeanEffectiveStress::deffective(const double *const s,
                                    double *const deff) const {

  for (int i = 0; i < 6; i++) {
    if (i < 3)
      deff[i] = 1. / 3.;
    else
      deff[i] = 0;
  }
  return 0;
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

int HuddlestonEffectiveStress::effective(const double * const s, double & eff) const
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
    return 0;
  }


  eff = se * exp(b_*(I1v / ss - 1.0));

  return 0;
}

int HuddlestonEffectiveStress::deffective(const double * const s, double * const deff) const
{
  // Useful common factors
  double sd[6];
  std::copy(s, s+6, sd);
  dev_vec(sd);
  std::fill(deff, deff+6, 0.0);
  
  if (norm2_vec(s, 6) == 0.0) 
    return 0;

  double I1v = I1(s);
  double I2v = I2(s);
  double I2pv = I2(sd);
  double se = sqrt(-3.0 * I2pv);
  double ss = sqrt(-3.0 * I2pv + I2v);

  // check if ss = 0.0
  if (ss == 0.0)
    return 0;
 
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

  return 0;
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

int MaxPrincipalEffectiveStress::effective(const double * const s, double & eff) const
{
  double vals[3];

  int ier = eigenvalues_sym(s, vals);
  eff = vals[2];

  if (eff < 0.0) eff = 0.0;

  return ier;
}

int MaxPrincipalEffectiveStress::deffective(const double * const s, double * const deff) const
{
  double vectors[9];
  double values[3];
  int ier = eigenvalues_sym(s, values);

  if (values[2] < 0.0) {
    std::fill(deff, deff+6, 0.0);
    return 0.0;
  }

  ier = eigenvectors_sym(s, vectors);
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

  return ier;
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

int MaxSeveralEffectiveStress::effective(const double * const s, double & eff) const
{
  size_t ind;
  select_(s, ind, eff);
  return 0;
}

int MaxSeveralEffectiveStress::deffective(const double * const s, double * const deff) const
{
  size_t ind;
  double eff;
  select_(s, ind, eff);

  measures_[ind]->deffective(s, deff);

  return 0;
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

int SumSeveralEffectiveStress::effective(const double * const s, double & eff) const
{
  double val;
  eff = 0.0;
  for (size_t i = 0; i < measures_.size(); i++) {
    measures_[i]->effective(s, val);
    eff += weights_[i] * val;
  }
  return 0;
}

int SumSeveralEffectiveStress::deffective(const double * const s, double * const deff) const
{
  double val[6];
  std::fill(deff, deff+6, 0.0);
  for (size_t i = 0; i < measures_.size(); i++) {
    measures_[i]->deffective(s, val);
    for (size_t j = 0; j < 6; j++) {
      deff[j] += weights_[i] * val[j];
    }
  }
  return 0;
}

ModularCreepDamageModel_sd::ModularCreepDamageModel_sd(ParameterSet & params) :
      NEMLScalarDamagedModel_sd(params),
      A_(params.get_object_parameter<Interpolate>("A")), 
      xi_(params.get_object_parameter<Interpolate>("xi")),
      phi_(params.get_object_parameter<Interpolate>("phi")), 
      estress_(params.get_object_parameter<EffectiveStress>("estress"))
{

}

std::string ModularCreepDamageModel_sd::type()
{
  return "ModularCreepDamageModel_sd";
}

ParameterSet ModularCreepDamageModel_sd::parameters()
{
  ParameterSet pset(ModularCreepDamageModel_sd::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("xi");
  pset.add_parameter<NEMLObject>("phi");
  pset.add_parameter<NEMLObject>("estress");
  pset.add_parameter<NEMLObject>("base");

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

std::unique_ptr<NEMLObject> ModularCreepDamageModel_sd::initialize(ParameterSet & params)
{
  return neml::make_unique<ModularCreepDamageModel_sd>(params); 
}

int ModularCreepDamageModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se;
  estress_->effective(s_np1, se);
  double dt = t_np1 - t_n;
  
  *dd = d_n + pow(se / A, xi) * pow(1.0 - d_np1, xi-phi) * dt;

  return 0;
}

int ModularCreepDamageModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se;
  estress_->effective(s_np1, se);
  double dt = t_np1 - t_n;

  *dd = pow(se / A, xi) * (phi-xi) * pow(1.0 - d_np1, xi-phi - 1.0) * dt;

  return 0;
}

int ModularCreepDamageModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);

  return 0;
}

int ModularCreepDamageModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double xi = xi_->value(T_np1);
  double A = A_->value(T_np1);
  double phi = phi_->value(T_np1);

  double se;
  estress_->effective(s_np1, se);
  double dt = t_np1 - t_n;

  if (se == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  double scalar = pow(se / A, xi - 1.0) * xi / A * pow(1.0 - d_np1, xi-phi) * dt;
  
  estress_->deffective(s_np1, dd);
  for (int i=0; i<6; i++) dd[i] *= scalar;
  
  return 0;
}

LarsonMillerCreepDamageModel_sd::LarsonMillerCreepDamageModel_sd(ParameterSet &
                                                                 params) :
      NEMLScalarDamagedModel_sd(params),
      lmr_(params.get_object_parameter<LarsonMillerRelation>("lmr")), 
      estress_(params.get_object_parameter<EffectiveStress>("estress"))
{
}

std::string LarsonMillerCreepDamageModel_sd::type()
{
  return "LarsonMillerCreepDamageModel_sd";
}

ParameterSet LarsonMillerCreepDamageModel_sd::parameters()
{
  ParameterSet pset(LarsonMillerCreepDamageModel_sd::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("lmr");
  pset.add_parameter<NEMLObject>("estress");
  pset.add_parameter<NEMLObject>("base");

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

std::unique_ptr<NEMLObject> LarsonMillerCreepDamageModel_sd::initialize(ParameterSet & params)
{
  return neml::make_unique<LarsonMillerCreepDamageModel_sd>(params); 
}

int LarsonMillerCreepDamageModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double se;
  estress_->effective(s_np1, se);
  double dt = t_np1 - t_n;

  // The usual problem
  if (se == 0.0) {
    *dd = d_n;
    return 0;
  }

  double tR;
  int ier = lmr_->tR(se * (1-d_np1), T_np1, tR);
  if (ier != 0) return ier;
  
  *dd = d_n + dt / tR;

  return 0;
}

int LarsonMillerCreepDamageModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double se;
  estress_->effective(s_np1, se);

  // The usual problem
  if (se == 0.0) {
    *dd = 0.0;
    return 0;
  }

  double dt = t_np1 - t_n;

  double tR;
  int ier = lmr_->tR(se * (1-d_np1), T_np1, tR);
  if (ier != 0) return ier;

  double dtR;
  ier = lmr_->dtR_ds(se * (1-d_np1), T_np1, dtR);
  if (ier != 0) return ier;

  *dd = (dt * dtR * se) / (tR * tR);
  
  return 0;
}

int LarsonMillerCreepDamageModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  std::fill(dd, dd+6, 0.0);

  return 0;
}

int LarsonMillerCreepDamageModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double se;
  estress_->effective(s_np1, se);
  double dt = t_np1 - t_n;

  // The usual problem
  if (se == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  double tR;
  int ier = lmr_->tR(se * (1-d_np1), T_np1, tR);
  if (ier != 0) return ier;

  double dtR;
  ier = lmr_->dtR_ds(se * (1-d_np1), T_np1, dtR);
  if (ier != 0) return ier;

  double dse = -dtR * dt * (1-d_np1) / (tR * tR);

  estress_->deffective(s_np1, dd);
  for (int i=0; i<6; i++) dd[i] *= dse;

  return 0;
}

NEMLStandardScalarDamagedModel_sd::NEMLStandardScalarDamagedModel_sd(ParameterSet
                                                                     & params) :
      NEMLScalarDamagedModel_sd(params) 
{

}

int NEMLStandardScalarDamagedModel_sd::damage(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  int ier = f(s_np1, d_np1, T_np1, fval);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);
  *dd = d_n + fval * deps;

  return ier;
}

int NEMLStandardScalarDamagedModel_sd::ddamage_dd(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double df;
  int ier = df_dd(s_np1, d_np1, T_np1, df);
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  *dd = df * deps;

  return ier;
}

int NEMLStandardScalarDamagedModel_sd::ddamage_de(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  int ier = f(s_np1, d_np1, T_np1, fval);
  if (ier != SUCCESS) return ier;
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  if (deps == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }
  
  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] - s_n[i];
    de[i] = e_np1[i] - e_n[i];
  }

  double S[36];
  ier = elastic_->S(T_np1, S);
  if (ier != SUCCESS) return ier;

  double dee[36];
  ier = mat_vec(S, 6, ds, 6, dee);
  if (ier != SUCCESS) return ier;

  for (int i=0; i<6; i++) {
    dd[i] = (2.0 * fval) / (3.0 * deps) * (de[i] - dee[i]); 
  }

  return 0;
}

int NEMLStandardScalarDamagedModel_sd::ddamage_ds(
    double d_np1, double d_n, 
    const double * const e_np1, const double * const e_n,
    const double * const s_np1, const double * const s_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const dd) const
{
  double fval;
  int ier = f(s_np1, d_np1, T_np1, fval);
  if (ier != SUCCESS) return ier;
  double deps = dep(s_np1, s_n, e_np1, e_n, T_np1);

  if (deps == 0.0) {
    std::fill(dd, dd+6, 0.0);
    return 0;
  }

  double ds[6];
  double de[6];
  for (int i=0; i<6; i++) {
    ds[i] = s_np1[i] - s_n[i];
    de[i] = e_np1[i] - e_n[i];
  }

  double S[36];
  ier = elastic_->S(T_np1, S);
  if (ier != SUCCESS) return ier;

  double dee[36];
  ier = mat_vec(S, 6, ds, 6, dee);
  if (ier != SUCCESS) return ier;

  double v1[6];
  for (int i=0; i<6; i++) {
    v1[i] = (2.0 * fval) / (3.0 * deps) * (dee[i] - de[i]); 
  }

  ier = mat_vec(S, 6, v1, 6, dd);
  if (ier != SUCCESS) return ier;

  double dds[6];
  ier = df_ds(s_np1, d_np1, T_np1, dds);
  if (ier != SUCCESS) return ier;

  for (int i=0; i<6; i++) {
    dd[i] = dds[i] * deps + dd[i];
  }

  return 0;
}

double NEMLStandardScalarDamagedModel_sd::dep(
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

NEMLWorkDamagedModel_sd::NEMLWorkDamagedModel_sd(ParameterSet & params) :
      NEMLScalarDamagedModel_sd(params), 
      Wcrit_(params.get_object_parameter<Interpolate>("Wcrit")),
      n_(params.get_parameter<double>("n")),
      eps_(params.get_parameter<double>("eps"))
{
}

std::string NEMLWorkDamagedModel_sd::type()
{
  return "NEMLWorkDamagedModel_sd";
}

ParameterSet NEMLWorkDamagedModel_sd::parameters()
{
  ParameterSet pset(NEMLWorkDamagedModel_sd::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("Wcrit");
  pset.add_parameter<double>("n");
  pset.add_parameter<NEMLObject>("base");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          make_constant(0.0));
  pset.add_optional_parameter<double>("rtol", 1.0e-14);
  pset.add_optional_parameter<double>("atol", 1.0e-8);
  pset.add_optional_parameter<int>("miter", 50);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<bool>("linesearch", false);

  pset.add_optional_parameter<bool>("truesdell", true);
  pset.add_optional_parameter<bool>("ekill", false);
  pset.add_optional_parameter<double>("dkill", 0.5);
  pset.add_optional_parameter<double>("sfact", 100000.0);

  pset.add_optional_parameter<double>("eps", 1e-30);

  return pset;
}

std::unique_ptr<NEMLObject> NEMLWorkDamagedModel_sd::initialize(ParameterSet & params)
{
  return neml::make_unique<NEMLWorkDamagedModel_sd>(params); 
}

int NEMLWorkDamagedModel_sd::damage(double d_np1, double d_n,
                   const double * const e_np1, const double * const e_n,
                   const double * const s_np1, const double * const s_n,
                   double T_np1, double T_n,
                   double t_np1, double t_n,
                   double * const dd) const
{
  if (d_np1 == 0.0) {
    *dd = d_n;
    return 0;
  }

  double wrate = workrate(e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n,
                          std::fabs(d_np1), d_n);
  if (wrate == 0.0) {
    *dd = d_n;
    return 0;
  }

  double dt = t_np1 - t_n;

  double val = Wcrit_->value(wrate);

  *dd = d_n + n_ * std::pow(std::fabs(d_np1), (n_-1.0)/n_) * 
      wrate * dt / val;
  
  return 0;
}

int NEMLWorkDamagedModel_sd::ddamage_dd(double d_np1, double d_n,
                   const double * const e_np1, const double * const e_n,
                   const double * const s_np1, const double * const s_n,
                   double T_np1, double T_n,
                   double t_np1, double t_n,
                   double * const dd) const
{
  if (d_np1 == 0.0) {
    *dd = d_n;
    return 0;
  }

  double wrate = workrate(e_np1, e_n, s_np1, s_n, T_np1, T_n, t_np1, t_n,
                          std::fabs(d_np1), d_n);
  if (wrate == 0.0) {
    *dd = d_n;
    return 0;
  }

  double val = Wcrit_->value(wrate);
  double deriv = Wcrit_->derivative(wrate);
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

  return 0;
}

int NEMLWorkDamagedModel_sd::ddamage_de(double d_np1, double d_n,
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
    return 0;
  }

  double val = Wcrit_->value(wrate);
  double dval = Wcrit_->derivative(wrate);

  double fact = n_ * std::pow(std::fabs(d_np1), (n_-1.0)/n_) / val * 
      (1.0 - wrate / val * dval) * (1.0 - std::fabs(d_np1));

  for (size_t i = 0; i < 6; i++) {
    dd[i] = fact * s_np1[i];
  }

  return 0;
}

int NEMLWorkDamagedModel_sd::ddamage_ds(double d_np1, double d_n,
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
    return 0;
  }

  double val = Wcrit_->value(wrate);
  double dval = Wcrit_->derivative(wrate);

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

  return 0;
}

double NEMLWorkDamagedModel_sd::workrate(
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

NEMLPowerLawDamagedModel_sd::NEMLPowerLawDamagedModel_sd(ParameterSet & params) :
      NEMLStandardScalarDamagedModel_sd(params), 
      A_(params.get_object_parameter<Interpolate>("A")),
      a_(params.get_object_parameter<Interpolate>("a"))
{

}

std::string NEMLPowerLawDamagedModel_sd::type()
{
  return "NEMLPowerLawDamagedModel_sd";
}

ParameterSet NEMLPowerLawDamagedModel_sd::parameters()
{
  ParameterSet pset(NEMLPowerLawDamagedModel_sd::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("a");
  pset.add_parameter<NEMLObject>("base");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          make_constant(0.0));
  pset.add_optional_parameter<double>("rtol", 1.0e-10);
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

std::unique_ptr<NEMLObject> NEMLPowerLawDamagedModel_sd::initialize(ParameterSet & params)
{
  return neml::make_unique<NEMLPowerLawDamagedModel_sd>(params); 
}

int NEMLPowerLawDamagedModel_sd::f(const double * const s_np1, double d_np1,
                                double T_np1, double & f) const
{
  double sev = se(s_np1);
  double A = A_->value(T_np1);
  double a = a_->value(T_np1);

  f = A * pow(sev, a);

  return 0;
}

int NEMLPowerLawDamagedModel_sd::df_ds(const double * const s_np1, double d_np1, double T_np1,
                                 double * const df) const
{
  double sev = se(s_np1);
  double A = A_->value(T_np1);
  double a = a_->value(T_np1);

  if (sev == 0.0) {
    std::fill(df, df+6, 0.0);
    return 0;
  }

  std::copy(s_np1, s_np1+6, df);
  double sm = (s_np1[0] + s_np1[1] + s_np1[2]) / 3.0;
  for (int i=0; i<3; i++) df[i] -= sm;
  for (int i=0; i<6; i++) df[i] *= (3.0 * A * a / 2.0 * pow(sev, a - 2.0));

  return 0;
}

int NEMLPowerLawDamagedModel_sd::df_dd(const double * const s_np1, double d_np1, double T_np1,
                                 double & df) const
{
  df = 0.0;

  return 0;
}

double NEMLPowerLawDamagedModel_sd::se(const double * const s) const
{
  return sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) + 
               pow(s[2] - s[0], 2.0) + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + 
                                              pow(s[5], 2.0))) / 2.0);
}

NEMLExponentialWorkDamagedModel_sd::NEMLExponentialWorkDamagedModel_sd(ParameterSet
                                                                       & params) :
      NEMLStandardScalarDamagedModel_sd(params), 
      W0_(params.get_object_parameter<Interpolate>("W0")),
      k0_(params.get_object_parameter<Interpolate>("k0")),
      af_(params.get_object_parameter<Interpolate>("af"))
{

}

std::string NEMLExponentialWorkDamagedModel_sd::type()
{
  return "NEMLExponentialWorkDamagedModel_sd";
}

ParameterSet NEMLExponentialWorkDamagedModel_sd::parameters()
{
  ParameterSet pset(NEMLExponentialWorkDamagedModel_sd::type());

  pset.add_parameter<NEMLObject>("elastic");
  pset.add_parameter<NEMLObject>("W0");
  pset.add_parameter<NEMLObject>("k0");
  pset.add_parameter<NEMLObject>("af");
  pset.add_parameter<NEMLObject>("base");

  pset.add_optional_parameter<NEMLObject>("alpha",
                                          make_constant(0.0));
  pset.add_optional_parameter<double>("rtol", 1.0e-10);
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

std::unique_ptr<NEMLObject> NEMLExponentialWorkDamagedModel_sd::initialize(ParameterSet & params)
{
  return neml::make_unique<NEMLExponentialWorkDamagedModel_sd>(params); 
}

int NEMLExponentialWorkDamagedModel_sd::f(const double * const s_np1, double d_np1,
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

  return 0;
}

int NEMLExponentialWorkDamagedModel_sd::df_ds(const double * const s_np1, double d_np1, double T_np1,
                                 double * const df) const
{
  double sev = se(s_np1);
  double W0 = W0_->value(T_np1);
  double k0 = k0_->value(T_np1);
  double af = af_->value(T_np1);

  if (sev == 0.0) {
    std::fill(df, df+6, 0.0);
    return 0;
  }

  if ((d_np1 + k0) < 0.0) {
    std::fill(df, df+6, 0.0);
    return 0;
  }

  std::copy(s_np1, s_np1+6, df);
  double sm = (s_np1[0] + s_np1[1] + s_np1[2]) / 3.0;
  for (int i=0; i<3; i++) df[i] -= sm;
  for (int i=0; i<6; i++) df[i] *= (3.0 * pow(d_np1 + k0, af) / (2.0 * sev * W0) );

  return 0;
}

int NEMLExponentialWorkDamagedModel_sd::df_dd(const double * const s_np1, double d_np1, double T_np1,
                                 double & df) const
{
  double sev = se(s_np1);
  double W0 = W0_->value(T_np1);
  double k0 = k0_->value(T_np1);
  double af = af_->value(T_np1);

  if ((d_np1 + k0) < 0.0) {
    df = 0.0;
    return 0;
  }

  df = af * pow(d_np1 + k0, af - 1.0) * sev / W0;

  return 0;
}

double NEMLExponentialWorkDamagedModel_sd::se(const double * const s) const
{
  return sqrt((pow(s[0]-s[1], 2.0) + pow(s[1] - s[2], 2.0) + 
               pow(s[2] - s[0], 2.0) + 3.0 * (pow(s[3], 2.0) + pow(s[4], 2.0) + 
                                              pow(s[5], 2.0))) / 2.0);
}

}
