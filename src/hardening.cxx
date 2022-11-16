#include "hardening.h"

#include "math/nemlmath.h"
#include "nemlerror.h"

#include <cmath>
#include <algorithm>
#include <iostream>

namespace neml {

HardeningRule::HardeningRule(ParameterSet & params) :
    HistoryNEMLObject(params)
{

}

IsotropicHardeningRule::IsotropicHardeningRule(ParameterSet & params) :
    HardeningRule(params)
{

}

void IsotropicHardeningRule::populate_hist(History & h) const
{
  h.add<double>(prefix("alpha"));
}

void IsotropicHardeningRule::init_hist(History & h) const
{
  h.get<double>(prefix("alpha")) = 0.0;
}

// Implementation of linear hardening
LinearIsotropicHardeningRule::LinearIsotropicHardeningRule(ParameterSet & params) :
    IsotropicHardeningRule(params),
    s0_(params.get_object_parameter<Interpolate>("s0")),
    K_(params.get_object_parameter<Interpolate>("K"))
{

}

std::string LinearIsotropicHardeningRule::type()
{
  return "LinearIsotropicHardeningRule";
}

ParameterSet LinearIsotropicHardeningRule::parameters()
{
  ParameterSet pset(LinearIsotropicHardeningRule::type());

  pset.add_parameter<NEMLObject>("s0");
  pset.add_parameter<NEMLObject>("K");

  return pset;
}

std::unique_ptr<NEMLObject> LinearIsotropicHardeningRule::initialize(ParameterSet & params)
{
  return neml::make_unique<LinearIsotropicHardeningRule>(params); 
}

void LinearIsotropicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  qv[0] = -s0_->value(T) - K_->value(T) * alpha[0];

}

void LinearIsotropicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  dqv[0] = -K_->value(T);

}

double LinearIsotropicHardeningRule::s0(double T) const
{
  return s0_->value(T);
}

double LinearIsotropicHardeningRule::K(double T) const
{
  return K_->value(T);
}


InterpolatedIsotropicHardeningRule::InterpolatedIsotropicHardeningRule(ParameterSet & params) :
    IsotropicHardeningRule(params),
    flow_(params.get_object_parameter<Interpolate>("flow"))
{

}

std::string InterpolatedIsotropicHardeningRule::type()
{
  return "InterpolatedIsotropicHardeningRule";
}

ParameterSet InterpolatedIsotropicHardeningRule::parameters()
{
  ParameterSet pset(InterpolatedIsotropicHardeningRule::type());

  pset.add_parameter<NEMLObject>("flow");

  return pset;
}

std::unique_ptr<NEMLObject> InterpolatedIsotropicHardeningRule::initialize(ParameterSet & params)
{
  return neml::make_unique<InterpolatedIsotropicHardeningRule>(
      params
      ); 
}

void InterpolatedIsotropicHardeningRule::q(const double * const alpha,
                                          double T, double * const qv) const
{
  qv[0] = -flow_->value(alpha[0]);
}

void InterpolatedIsotropicHardeningRule::dq_da(const double * const alpha,
                                              double T,
                                              double * const dqv) const
{
  dqv[0] = -flow_->derivative(alpha[0]);
}

// Implementation of voce hardening
VoceIsotropicHardeningRule::VoceIsotropicHardeningRule(
    ParameterSet & params) :
      IsotropicHardeningRule(params),
      s0_(params.get_object_parameter<Interpolate>("s0")), 
      R_(params.get_object_parameter<Interpolate>("R")),
      d_(params.get_object_parameter<Interpolate>("d"))
{

}

std::string VoceIsotropicHardeningRule::type()
{
  return "VoceIsotropicHardeningRule";
}

ParameterSet VoceIsotropicHardeningRule::parameters()
{
  ParameterSet pset(VoceIsotropicHardeningRule::type());

  pset.add_parameter<NEMLObject>("s0");
  pset.add_parameter<NEMLObject>("R");
  pset.add_parameter<NEMLObject>("d");

  return pset;
}

std::unique_ptr<NEMLObject> VoceIsotropicHardeningRule::initialize(ParameterSet & params)
{
  return neml::make_unique<VoceIsotropicHardeningRule>(params); 
}

void VoceIsotropicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  qv[0] = -s0_->value(T) - R_->value(T) * (1.0 - exp(-d_->value(T) * alpha[0]));

}

void VoceIsotropicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  dqv[0] = -d_->value(T) * R_->value(T) * exp(-d_->value(T) * alpha[0]);

}

double VoceIsotropicHardeningRule::s0(double T) const
{
  return s0_->value(T);
}

double VoceIsotropicHardeningRule::R(double T) const
{
  return R_->value(T);
}

double VoceIsotropicHardeningRule::d(double T) const
{
  return d_->value(T);
}

// Implementation of voce hardening
PowerLawIsotropicHardeningRule::PowerLawIsotropicHardeningRule(ParameterSet &
                                                               params) :
    IsotropicHardeningRule(params),
    s0_(params.get_object_parameter<Interpolate>("s0")),
    A_(params.get_object_parameter<Interpolate>("A")),
    n_(params.get_object_parameter<Interpolate>("n"))
{

}

std::string PowerLawIsotropicHardeningRule::type()
{
  return "PowerLawIsotropicHardeningRule";
}

ParameterSet PowerLawIsotropicHardeningRule::parameters()
{
  ParameterSet pset(PowerLawIsotropicHardeningRule::type());

  pset.add_parameter<NEMLObject>("s0");
  pset.add_parameter<NEMLObject>("A");
  pset.add_parameter<NEMLObject>("n");

  return pset;
}

std::unique_ptr<NEMLObject> PowerLawIsotropicHardeningRule::initialize(ParameterSet & params)
{
  return neml::make_unique<PowerLawIsotropicHardeningRule>(params); 
}

void PowerLawIsotropicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  qv[0] = -s0_->value(T) - A_->value(T) * pow(alpha[0], n_->value(T));

}

void PowerLawIsotropicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  if (alpha[0] == 0.0) {
    dqv[0] = -1e15; // In actuality - infinity
  }
  else {
    dqv[0] = -A_->value(T) * n_->value(T) * pow(alpha[0], n_->value(T) - 1);
  }

}

// Implementation of combined isotropic class
CombinedIsotropicHardeningRule::CombinedIsotropicHardeningRule(ParameterSet &
                                                               params): 
    IsotropicHardeningRule(params),    
    rules_(params.get_object_parameter_vector<IsotropicHardeningRule>("rules"))
{

}

std::string CombinedIsotropicHardeningRule::type()
{
  return "CombinedIsotropicHardeningRule";
}

ParameterSet CombinedIsotropicHardeningRule::parameters()
{
  ParameterSet pset(CombinedIsotropicHardeningRule::type());

  pset.add_parameter<std::vector<NEMLObject>>("rules");

  return pset;
}

std::unique_ptr<NEMLObject> CombinedIsotropicHardeningRule::initialize(ParameterSet & params)
{
  return neml::make_unique<CombinedIsotropicHardeningRule>(
      params
      ); 
}

void CombinedIsotropicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  double qi;
  qv[0] = 0.0;
  for (auto it = rules_.begin(); it != rules_.end(); ++it) {
    (*it)->q(alpha, T, &qi);
    qv[0] += qi;
  }
}

void CombinedIsotropicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  double dqi;
  dqv[0] = 0.0;
  for (auto it = rules_.begin(); it != rules_.end(); ++it) {
    (*it)->dq_da(alpha, T, &dqi);
    dqv[0] += dqi;
  }
}

size_t CombinedIsotropicHardeningRule::nrules() const 
{
  return rules_.size();
}

KinematicHardeningRule::KinematicHardeningRule(ParameterSet & params) :
    HardeningRule(params)
{

}

// Implementation of kinematic base class
void KinematicHardeningRule::populate_hist(History & hist) const
{
  hist.add<Symmetric>(prefix("backstress"));
}

void KinematicHardeningRule::init_hist(History & hist) const
{
  hist.get<Symmetric>(prefix("backstress")) = Symmetric::zero();
}

// Implementation of linear kinematic hardening
LinearKinematicHardeningRule::LinearKinematicHardeningRule(ParameterSet & params) :
    KinematicHardeningRule(params),
    H_(params.get_object_parameter<Interpolate>("H"))
{

}

std::string LinearKinematicHardeningRule::type()
{
  return "LinearKinematicHardeningRule";
}

ParameterSet LinearKinematicHardeningRule::parameters()
{
  ParameterSet pset(LinearKinematicHardeningRule::type());

  pset.add_parameter<NEMLObject>("H");

  return pset;
}

std::unique_ptr<NEMLObject> LinearKinematicHardeningRule::initialize(ParameterSet & params)
{
  return neml::make_unique<LinearKinematicHardeningRule>(params); 
}

void LinearKinematicHardeningRule::q(const double * const alpha, 
                                    double T, double * const qv) const
{
  for (size_t i=0; i<6; i++) {
    qv[i] = -H_->value(T) * alpha[i];
  }

}

void LinearKinematicHardeningRule::dq_da(const double * const alpha, 
                                    double T, double * const dqv) const
{
  std::fill(dqv, dqv+36, 0.0);
  for (size_t i=0; i<6; i++) {
    dqv[CINDEX(i,i,6)] = -H_->value(T);
  }

}

double LinearKinematicHardeningRule::H(double T) const
{
  return H_->value(T);
}

CombinedHardeningRule::CombinedHardeningRule(ParameterSet & params) :
    HardeningRule(params),
    iso_(params.get_object_parameter<IsotropicHardeningRule>("iso")), 
    kin_(params.get_object_parameter<KinematicHardeningRule>("kin"))
{

}

std::string CombinedHardeningRule::type()
{
  return "CombinedHardeningRule";
}

ParameterSet CombinedHardeningRule::parameters()
{
  ParameterSet pset(CombinedHardeningRule::type());

  pset.add_parameter<NEMLObject>("iso");
  pset.add_parameter<NEMLObject>("kin");

  return pset;
}

std::unique_ptr<NEMLObject> CombinedHardeningRule::initialize(ParameterSet & params)
{
  return neml::make_unique<CombinedHardeningRule>(params); 
}

void CombinedHardeningRule::populate_hist(History & hist) const
{
  iso_->populate_hist(hist);
  kin_->populate_hist(hist);
}

void CombinedHardeningRule::init_hist(History & hist) const
{
  iso_->init_hist(hist);
  kin_->init_hist(hist);
}

void CombinedHardeningRule::q(const double * const alpha, double T, 
                             double * const qv) const
{
  iso_->q(alpha, T, qv);
  return kin_->q(&alpha[iso_->nhist()], T, &qv[iso_->nhist()]);
}

void CombinedHardeningRule::dq_da(const double * const alpha, double T, 
                                 double * const dqv) const
{
  // Annoying this doesn't work nicely...
  std::vector<double> idv(iso_->nhist() * iso_->nhist());
  double * id = &idv[0];
  iso_->dq_da(alpha, T, id);
  
  std::vector<double> kdv(kin_->nhist() * kin_->nhist());
  double * kd = &kdv[0];
  kin_->dq_da(&alpha[iso_->nhist()], T, kd);

  std::fill(dqv, dqv + nhist()*nhist(), 0.0);

  for (size_t i=0; i<iso_->nhist(); i++) {
    for (size_t j=0; j<iso_->nhist(); j++) {
      dqv[CINDEX(i,j,nhist())] = id[CINDEX(i,j,iso_->nhist())];
    }
  }
  
  size_t os = iso_->nhist();
  for (size_t i=0; i<kin_->nhist(); i++) {
    for (size_t j=0; j<kin_->nhist(); j++) {
      dqv[CINDEX((i+os),(j+os),nhist())] = kd[CINDEX(i,j,kin_->nhist())];
    }
  }

}

NonAssociativeHardening::NonAssociativeHardening(ParameterSet & params) :
    HistoryNEMLObject(params)
{

}

// Provide zeros for these
void NonAssociativeHardening::h_time(const double * const s, 
                                    const double * const alpha, double T,
                                    double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
}

void NonAssociativeHardening::dh_ds_time(const double * const s,
                                        const double * const alpha, double T,
                                        double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*6, 0.0);
}

void NonAssociativeHardening::dh_da_time(const double * const s,
                                        const double * const alpha, double T,
                                        double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);
}

// Provide zeros for these
void NonAssociativeHardening::h_temp(const double * const s,
                                    const double * const alpha, double T,
                                    double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
}

void NonAssociativeHardening::dh_ds_temp(const double * const s,
                                        const double * const alpha, double T,
                                        double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*6, 0.0);
}

void NonAssociativeHardening::dh_da_temp(const double * const s,
                                        const double * const alpha, double T,
                                        double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);
}

// Begin non-associative hardening rules
//
// Gamma functions for Chaboche
GammaModel::GammaModel(ParameterSet & params) :
    NEMLObject(params)
{

}
//
// Constant
//
ConstantGamma::ConstantGamma(ParameterSet & params) :
    GammaModel(params),
    g_(params.get_object_parameter<Interpolate>("g"))
{

}

std::string ConstantGamma::type()
{
  return "ConstantGamma";
}

ParameterSet ConstantGamma::parameters()
{
  ParameterSet pset(ConstantGamma::type());

  pset.add_parameter<NEMLObject>("g");

  return pset;
}

std::unique_ptr<NEMLObject> ConstantGamma::initialize(ParameterSet & params)
{
  return neml::make_unique<ConstantGamma>(
      params
      ); 
}

double ConstantGamma::gamma(double ep, double T) const {
  return g_->value(T);
}

double ConstantGamma::dgamma(double ep, double T) const {
  return 0.0;
}

double ConstantGamma::g(double T) const {
  return g_->value(T);
}

//
// Saturating
//
SatGamma::SatGamma(ParameterSet & params) :
    GammaModel(params),
    gs_(params.get_object_parameter<Interpolate>("gs")), 
    g0_(params.get_object_parameter<Interpolate>("g0")), 
    beta_(params.get_object_parameter<Interpolate>("beta"))
{

}

std::string SatGamma::type()
{
  return "SatGamma";
}

ParameterSet SatGamma::parameters()
{
  ParameterSet pset(SatGamma::type());

  pset.add_parameter<NEMLObject>("gs");
  pset.add_parameter<NEMLObject>("g0");
  pset.add_parameter<NEMLObject>("beta");

  return pset;
}

std::unique_ptr<NEMLObject> SatGamma::initialize(ParameterSet & params)
{
  return neml::make_unique<SatGamma>(params); 
}

double SatGamma::gamma(double ep, double T) const {
  return gs_->value(T) + (g0_->value(T) - gs_->value(T)) * exp(-beta_->value(T) * ep);
}

double SatGamma::dgamma(double ep, double T) const {
  return beta_->value(T) * (gs_->value(T) - g0_->value(T)) * exp(-beta_->value(T) * ep);
}

double SatGamma::gs(double T) const {
  return gs_->value(T);
}

double SatGamma::g0(double T) const {
  return g0_->value(T);
}

double SatGamma::beta(double T) const {
  return beta_->value(T);
}

//
// Chaboche
//
Chaboche::Chaboche(ParameterSet & params) :
    NonAssociativeHardening(params),
    iso_(params.get_object_parameter<IsotropicHardeningRule>("iso")),
    c_(params.get_object_parameter_vector<Interpolate>("C")),
    n_(c_.size()), 
    gmodels_(params.get_object_parameter_vector<GammaModel>("gmodels")), 
    A_(params.get_object_parameter_vector<Interpolate>("A")), 
    a_( params.get_object_parameter_vector<Interpolate>("a")), 
    relax_(true), 
    noniso_(params.get_parameter<bool>("noniso"))
{

}

std::string Chaboche::type()
{
  return "Chaboche";
}

ParameterSet Chaboche::parameters()
{
  ParameterSet pset(Chaboche::type());

  pset.add_parameter<NEMLObject>("iso");
  pset.add_parameter<std::vector<NEMLObject>>("C");
  pset.add_parameter<std::vector<NEMLObject>>("gmodels");
  pset.add_parameter<std::vector<NEMLObject>>("A");
  pset.add_parameter<std::vector<NEMLObject>>("a");

  pset.add_optional_parameter<bool>("noniso", true);

  return pset;
}

std::unique_ptr<NEMLObject> Chaboche::initialize(ParameterSet & params)
{
  return neml::make_unique<Chaboche>(params); 
}

size_t Chaboche::ninter() const
{
  return 1 + 6;
}

void Chaboche::populate_hist(History & h) const
{
  h.add<double>(prefix("alpha"));
  for (size_t i = 0; i < n_; i++)
    h.add<Symmetric>(prefix("backstress_"+std::to_string(i)));
}

void Chaboche::init_hist(History & h) const
{
  h.get<double>(prefix("alpha")) = 0.0;
  for (size_t i = 0; i < n_; i++)
    h.get<Symmetric>(prefix("backstress_"+std::to_string(i))) = Symmetric::zero();
}

void Chaboche::q(const double * const alpha, double T, double * const qv) const
{
  iso_->q(alpha, T, qv);
  std::fill(qv+1, qv+7, 0.0);
  
  // Helps with unrolling
  size_t n = n_;

  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<6; j++) {
      qv[j+1] += alpha[1+i*6+j];
    }
  }
}

void Chaboche::dq_da(const double * const alpha, double T, double * const qv) const
{
  std::fill(qv, qv+(ninter()*nhist()), 0.0);
  iso_->dq_da(alpha, T, qv); // fills in (0,0)
  
  // Help unroll
  size_t n = n_;

  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<6; j++) {
      qv[CINDEX((j+1),(1+i*6+j), nhist())] = 1.0;
    }
  }
}

void Chaboche::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  hv[0] = sqrt(2.0/3.0); // Isotropic part

  double X[6], nv[6];
  backstress_(alpha, X);
  std::copy(s, s+6, nv);
  dev_vec(nv);
  add_vec(nv, X, 6, nv);
  normalize_vec(nv, 6);
  
  // Note the extra factor of sqrt(2.0/3.0) -- this is to make it equivalent
  // to Chaboche's original definition
  
  std::vector<double> c = eval_vector(c_, T);

  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      hv[1+i*6+j] = - 2.0 / 3.0 * c[i] * nv[j] - sqrt(2.0/3.0) * gmodels_[i]->gamma(alpha[0], T) * alpha[1+i*6+j];
    }
  }

}

void Chaboche::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv + nhist()*6, 0.0);

  std::vector<double> c = eval_vector(c_, T);

  double X[6];
  backstress_(alpha, X);

  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, X, 6, n);
  double nv = norm2_vec(n, 6);
  normalize_vec(n, 6);
  
  double nn[36];

  std::fill(nn, nn+36, 0.0);
  for (size_t i=0; i<6; i++) {
    nn[CINDEX(i,i,6)] += 1.0;
  }
  
  double iv[6];
  double jv[6];
  for (size_t i=0; i<3; i++) {
    iv[i] = 1.0 / 3.0;
    jv[i] = 1.0;
  }
  for (int i=3; i<6; i++) {
    iv[i] = 0.0;
    jv[i] = 0.0;
  }

  outer_update_minus(iv, 6, jv, 6, nn);

  outer_update_minus(n, 6, n, 6, nn);
  if (nv != 0.0) {
    for (size_t i=0; i<36; i++) {
      nn[i] /= nv;
    }
  }
  
  // Fill in...
  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      for (int k=0; k<6; k++) {
        dhv[CINDEX((1+i*6+j),(k),6)] = -2.0 / 3.0 * c[i] * nn[CINDEX(j,k,6)];
      }
    }
  }
  
}

void Chaboche::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  // Again, there should be no earthly reason the compiler can't do this
  size_t nh = nhist();

  std::fill(dhv, dhv + nh*nh, 0.0);

  std::vector<double> c = eval_vector(c_, T);

  double X[6];
  backstress_(alpha, X);

  double ss[36];
  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, X, 6, n);
  double nv = norm2_vec(n, 6);
  normalize_vec(n, 6);
  
  std::fill(ss, ss+36, 0.0);
  for (size_t i=0; i<6; i++) {
    ss[CINDEX(i,i,6)] += 1.0;
  }
  
  outer_update_minus(n, 6, n, 6, ss);
  if (nv != 0.0) {
    for (size_t i=0; i<36; i++) {
      ss[i] /= nv;
    }
  }
  
  // Fill in the gamma part
  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      dhv[CINDEX((1+i*6+j),(1+i*6+j),nh)] -= sqrt(2.0/3.0) * gmodels_[i]->gamma(alpha[0], T);
    }
  }

  // Fill in the ss part
  for (size_t bi=0; bi<n_; bi++) {
    for (size_t i=0; i<6; i++) {
      for (size_t bj=0; bj<n_; bj++) {
        for (size_t j=0; j<6; j++) {
          dhv[CINDEX((1+bi*6+i),(1+bj*6+j),nh)] -= 2.0 / 3.0 * c[bi]  * 
              ss[CINDEX(i,j,6)];
        }
      }
    }
  }

  // Fill in the alpha part
  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      dhv[CINDEX((1+i*6+j),0,nhist())] = -sqrt(2.0/3.0) * 
          gmodels_[i]->dgamma(alpha[0], T) * alpha[1+i*6+j];
    }
  }

}

void Chaboche::h_time(const double * const s, const double * const alpha, 
                     double T, double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);
 
  std::vector<double> A = eval_vector(A_, T);
  std::vector<double> a = eval_vector(a_, T);

  double Xi[6];
  double nXi;
  for (size_t i=0; i<n_; i++) {
    std::copy(&alpha[1+i*6], &alpha[1+(i+1)*6], Xi);
    nXi = norm2_vec(Xi, 6);
    for (size_t j=0; j<6; j++) {
      hv[1+i*6+j] = -A[i] * sqrt(3.0/2.0) * pow(nXi, a[i] - 1.0) *
          alpha[1+i*6+j];
    }
  }

}

void Chaboche::dh_ds_time(const double * const s, const double * const alpha, 
                         double T, double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*6, 0.0);

  // Also return if relax

}

void Chaboche::dh_da_time(const double * const s, const double * const alpha,
                         double T, double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);

  std::vector<double> A = eval_vector(A_, T);
  std::vector<double> a = eval_vector(a_, T);

  size_t nh = nhist();
  size_t n = n_;
  
  double XX[36];
  double Xi[6];
  double nXi;
  int ia,ib;
  double d;
  for (size_t i=0; i<n; i++) {
    std::copy(&alpha[1+i*6], &alpha[1+(i+1)*6], Xi);
    nXi = norm2_vec(Xi, 6);
    normalize_vec(Xi, 6);
    outer_vec(Xi, 6, Xi, 6, XX);
    for (size_t j=0; j<6; j++) {
      ia = 1 + i*6 + j;
      for (size_t k=0; k<6; k++) {
        ib = 1 + i*6 + k;
        if (j == k) {
          d = 1.0;
        }
        else {
          d = 0.0;
        }
        dhv[CINDEX(ia,ib,nh)] = -A[i] * sqrt(3.0/2.0) * pow(nXi, a[i]-1.0) * (
            d + (a[i] - 1.0) * XX[CINDEX(j,k,6)]);
      }
    }
  }

}

void Chaboche::h_temp(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  std::fill(hv, hv+nhist(), 0.0);

  std::vector<double> c = eval_vector(c_, T);
  std::vector<double> dc = eval_deriv_vector(c_, T);

  for (size_t i=0; i<n_; i++) {
    if (c[i] == 0.0) continue;
    for (size_t j=0; j<6; j++) {
      hv[1+i*6+j] = -sqrt(2.0/3.0) * dc[i] / c[i] * alpha[1+i*6+j];
    }
  }
}

void Chaboche::dh_ds_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*6, 0.0);
}

void Chaboche::dh_da_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv+nhist()*nhist(), 0.0);

  std::vector<double> c = eval_vector(c_, T);
  std::vector<double> dc = eval_deriv_vector(c_, T);

  for (size_t i=0; i<n_; i++) {
    if (c[i] == 0.0) continue;
    for (size_t j=0; j<6; j++) {
      int ci = 1 + i*6 + j;
      dhv[CINDEX(ci,ci,nhist())] = - sqrt(2.0/3.0) * dc[i] / c[i];
    }
  }

}

int Chaboche::n() const
{
  return n_;
}

std::vector<double> Chaboche::c(double T) const
{
  return eval_vector(c_, T);
}

void Chaboche::backstress_(const double * const alpha, double * const X) const
{
  std::fill(X, X+6, 0.0);
  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      X[j] += alpha[1+i*6+j];
    }
  }
}

//
// ChabocheVoceRecovery with fancy, hard-coded Voce isotropic hardening
//
ChabocheVoceRecovery::ChabocheVoceRecovery(ParameterSet & params) :
    NonAssociativeHardening(params),
    s0_(params.get_object_parameter<Interpolate>("s0")),
    theta0_(params.get_object_parameter<Interpolate>("theta0")),
    Rmax_(params.get_object_parameter<Interpolate>("Rmax")),
    Rmin_(params.get_object_parameter<Interpolate>("Rmin")),
    r1_(params.get_object_parameter<Interpolate>("r1")),
    r2_(params.get_object_parameter<Interpolate>("r2")),
    c_(params.get_object_parameter_vector<Interpolate>("C")),
    n_(c_.size()), 
    gmodels_(params.get_object_parameter_vector<GammaModel>("gmodels")), 
    A_(params.get_object_parameter_vector<Interpolate>("A")), 
    a_( params.get_object_parameter_vector<Interpolate>("a")), 
    noniso_(params.get_parameter<bool>("noniso"))
{
  cache_history_();
}

std::string ChabocheVoceRecovery::type()
{
  return "ChabocheVoceRecovery";
}

ParameterSet ChabocheVoceRecovery::parameters()
{
  ParameterSet pset(ChabocheVoceRecovery::type());

  pset.add_parameter<NEMLObject>("s0");
  pset.add_parameter<NEMLObject>("theta0");
  pset.add_parameter<NEMLObject>("Rmax");
  pset.add_parameter<NEMLObject>("Rmin");
  pset.add_parameter<NEMLObject>("r1");
  pset.add_parameter<NEMLObject>("r2");

  pset.add_parameter<std::vector<NEMLObject>>("C");
  pset.add_parameter<std::vector<NEMLObject>>("gmodels");
  pset.add_parameter<std::vector<NEMLObject>>("A");
  pset.add_parameter<std::vector<NEMLObject>>("a");

  pset.add_optional_parameter<bool>("noniso", true);

  return pset;
}

std::unique_ptr<NEMLObject> ChabocheVoceRecovery::initialize(ParameterSet & params)
{
  return neml::make_unique<ChabocheVoceRecovery>(params); 
}

size_t ChabocheVoceRecovery::ninter() const
{
  return 1 + 6;
}

void ChabocheVoceRecovery::populate_hist(History & h) const
{
  h.add<double>(prefix("alpha"));
  for (size_t i = 0; i < n_; i++)
    h.add<Symmetric>(prefix("backstress_"+std::to_string(i)));
}

void ChabocheVoceRecovery::init_hist(History & h) const
{
  h.get<double>(prefix("alpha")) = 0.0;
  for (size_t i = 0; i < n_; i++)
    h.get<Symmetric>(prefix("backstress_"+std::to_string(i))) = Symmetric::zero();
}

void ChabocheVoceRecovery::q(const double * const alpha, double T, double * const qv) const
{
  // Isotropic part
  qv[0] = -(s0_->value(T) + alpha[0]);

  std::fill(qv+1, qv+7, 0.0);
  
  // Helps with unrolling
  size_t n = n_;

  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<6; j++) {
      qv[j+1] += alpha[1+i*6+j];
    }
  }
}

void ChabocheVoceRecovery::dq_da(const double * const alpha, double T, double * const qv) const
{
  std::fill(qv, qv+(ninter()*nh()), 0.0);

  // Isotropic part
  qv[0] = -1.0;
  
  // Help unroll
  size_t n = n_;

  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<6; j++) {
      qv[CINDEX((j+1),(1+i*6+j), nh())] = 1.0;
    }
  }
}

void ChabocheVoceRecovery::h(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  // Isotropic
  hv[0] = theta0_->value(T) * (1 - alpha[0]/Rmax_->value(T)) *
      std::sqrt(2.0/3.0);

  double X[6], nv[6];
  backstress_(alpha, X);
  std::copy(s, s+6, nv);
  dev_vec(nv);
  add_vec(nv, X, 6, nv);
  normalize_vec(nv, 6);
  
  // Note the extra factor of sqrt(2.0/3.0) -- this is to make it equivalent
  // to Chaboche's original definition
  
  std::vector<double> c = eval_vector(c_, T);

  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      hv[1+i*6+j] = - 2.0 / 3.0 * c[i] * nv[j] - sqrt(2.0/3.0) * gmodels_[i]->gamma(alpha[0], T) * alpha[1+i*6+j];
    }
  }

}

void ChabocheVoceRecovery::dh_ds(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv + nh()*6, 0.0);

  std::vector<double> c = eval_vector(c_, T);

  double X[6];
  backstress_(alpha, X);

  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, X, 6, n);
  double nv = norm2_vec(n, 6);
  normalize_vec(n, 6);
  
  double nn[36];

  std::fill(nn, nn+36, 0.0);
  for (size_t i=0; i<6; i++) {
    nn[CINDEX(i,i,6)] += 1.0;
  }
  
  double iv[6];
  double jv[6];
  for (size_t i=0; i<3; i++) {
    iv[i] = 1.0 / 3.0;
    jv[i] = 1.0;
  }
  for (int i=3; i<6; i++) {
    iv[i] = 0.0;
    jv[i] = 0.0;
  }

  outer_update_minus(iv, 6, jv, 6, nn);

  outer_update_minus(n, 6, n, 6, nn);
  if (nv != 0.0) {
    for (size_t i=0; i<36; i++) {
      nn[i] /= nv;
    }
  }
  
  // Fill in...
  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      for (int k=0; k<6; k++) {
        dhv[CINDEX((1+i*6+j),(k),6)] = -2.0 / 3.0 * c[i] * nn[CINDEX(j,k,6)];
      }
    }
  }
  
}

void ChabocheVoceRecovery::dh_da(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  // Again, there should be no earthly reason the compiler can't do this
  size_t nhi = nh();

  std::fill(dhv, dhv + nhi*nhi, 0.0);

  // Isotropic contribution
  dhv[0] = -theta0_->value(T) / Rmax_->value(T) * std::sqrt(2.0/3.0);
  std::vector<double> c = eval_vector(c_, T);

  double X[6];
  backstress_(alpha, X);

  double ss[36];
  double n[6];
  std::copy(s, s+6, n);
  dev_vec(n);
  add_vec(n, X, 6, n);
  double nv = norm2_vec(n, 6);
  normalize_vec(n, 6);
  
  std::fill(ss, ss+36, 0.0);
  for (size_t i=0; i<6; i++) {
    ss[CINDEX(i,i,6)] += 1.0;
  }
  
  outer_update_minus(n, 6, n, 6, ss);
  if (nv != 0.0) {
    for (size_t i=0; i<36; i++) {
      ss[i] /= nv;
    }
  }
  
  // Fill in the gamma part
  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      dhv[CINDEX((1+i*6+j),(1+i*6+j),nhi)] -= sqrt(2.0/3.0) * gmodels_[i]->gamma(alpha[0], T);
    }
  }

  // Fill in the ss part
  for (size_t bi=0; bi<n_; bi++) {
    for (size_t i=0; i<6; i++) {
      for (size_t bj=0; bj<n_; bj++) {
        for (size_t j=0; j<6; j++) {
          dhv[CINDEX((1+bi*6+i),(1+bj*6+j),nhi)] -= 2.0 / 3.0 * c[bi]  * 
              ss[CINDEX(i,j,6)];
        }
      }
    }
  }

  // Fill in the alpha part
  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      dhv[CINDEX((1+i*6+j),0,nhi)] = -sqrt(2.0/3.0) * 
          gmodels_[i]->dgamma(alpha[0], T) * alpha[1+i*6+j];
    }
  }

}

void ChabocheVoceRecovery::h_time(const double * const s, const double * const alpha, 
                     double T, double * const hv) const
{
  std::fill(hv, hv+nh(), 0.0);

  // Isotropic recovery term 
  hv[0] = r1_->value(T) * (Rmin_->value(T) - alpha[0]
                           ) * std::pow(std::fabs(Rmin_->value(T) - alpha[0]),
                                        r2_->value(T) - 1.0);
 
  std::vector<double> A = eval_vector(A_, T);
  std::vector<double> a = eval_vector(a_, T);

  double Xi[6];
  double nXi;
  for (size_t i=0; i<n_; i++) {
    std::copy(&alpha[1+i*6], &alpha[1+(i+1)*6], Xi);
    nXi = norm2_vec(Xi, 6);
    for (size_t j=0; j<6; j++) {
      hv[1+i*6+j] = -A[i] * pow(sqrt(3.0/2.0) * nXi, a[i] - 1.0) *
          alpha[1+i*6+j];
    }
  }

}

void ChabocheVoceRecovery::dh_ds_time(const double * const s, const double * const alpha, 
                         double T, double * const dhv) const
{
  std::fill(dhv, dhv+nh()*6, 0.0);

  // Also return if relax

}

void ChabocheVoceRecovery::dh_da_time(const double * const s, const double * const alpha,
                         double T, double * const dhv) const
{
  std::fill(dhv, dhv+nh()*nh(), 0.0);

  dhv[0] = std::copysign(r1_->value(T) * r2_->value(T) *
                         std::pow(std::fabs(Rmin_->value(T) - alpha[0]), r2_->value(T) -
                                  1.0), Rmin_->value(T) - alpha[0]);

  std::vector<double> A = eval_vector(A_, T);
  std::vector<double> a = eval_vector(a_, T);

  size_t nhi = nh();
  size_t n = n_;
  
  double XX[36];
  double Xi[6];
  double nXi;
  int ia,ib;
  double d;
  for (size_t i=0; i<n; i++) {
    std::copy(&alpha[1+i*6], &alpha[1+(i+1)*6], Xi);
    nXi = norm2_vec(Xi, 6);
    normalize_vec(Xi, 6);
    outer_vec(Xi, 6, Xi, 6, XX);
    for (size_t j=0; j<6; j++) {
      ia = 1 + i*6 + j;
      for (size_t k=0; k<6; k++) {
        ib = 1 + i*6 + k;
        if (j == k) {
          d = 1.0;
        }
        else {
          d = 0.0;
        }
        dhv[CINDEX(ia,ib,nhi)] = -A[i] * pow( sqrt(3.0/2.0) * nXi, a[i]-1.0) * (
            d + (a[i] - 1.0) * XX[CINDEX(j,k,6)]);
      }
    }
  }

}

void ChabocheVoceRecovery::h_temp(const double * const s, const double * const alpha, double T,
              double * const hv) const
{
  std::fill(hv, hv+nh(), 0.0);

  std::vector<double> c = eval_vector(c_, T);
  std::vector<double> dc = eval_deriv_vector(c_, T);

  for (size_t i=0; i<n_; i++) {
    if (c[i] == 0.0) continue;
    for (size_t j=0; j<6; j++) {
      hv[1+i*6+j] = -sqrt(2.0/3.0) * dc[i] / c[i] * alpha[1+i*6+j];
    }
  }
}

void ChabocheVoceRecovery::dh_ds_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv+nh()*6, 0.0);
}

void ChabocheVoceRecovery::dh_da_temp(const double * const s, const double * const alpha, double T,
              double * const dhv) const
{
  std::fill(dhv, dhv+nh()*nh(), 0.0);

  std::vector<double> c = eval_vector(c_, T);
  std::vector<double> dc = eval_deriv_vector(c_, T);

  for (size_t i=0; i<n_; i++) {
    if (c[i] == 0.0) continue;
    for (size_t j=0; j<6; j++) {
      size_t ci = 1 + i*6 + j;
      dhv[CINDEX(ci,ci,nh())] = - sqrt(2.0/3.0) * dc[i] / c[i];
    }
  }
}

int ChabocheVoceRecovery::n() const
{
  return n_;
}

void ChabocheVoceRecovery::backstress_(const double * const alpha, double * const X) const
{
  std::fill(X, X+6, 0.0);
  for (size_t i=0; i<n_; i++) {
    for (size_t j=0; j<6; j++) {
      X[j] += alpha[1+i*6+j];
    }
  }
}

} // namespace neml
