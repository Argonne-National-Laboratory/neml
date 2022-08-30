#include "cp/slipharden.h"

#include <stdexcept>

namespace neml {

SlipHardening::SlipHardening(ParameterSet & params): 
    HistoryNEMLObject(params)
{

}

History SlipHardening::cache(CacheType type) const
{
  switch (type) {
    case CacheType::BLANK: return *blank_;
      break;
    case CacheType::DOUBLE: return *double_;
      break;
    default: throw std::runtime_error("Invalid cached type");
      break;
  }
}

void SlipHardening::init_cache_()
{
  blank_ = make_unique<History>();
  populate_hist(*blank_);
  blank_->zero();

  double_ = make_unique<History>((*blank_).derivative<double>());
}

History SlipHardening::d_hist_to_tau_ext(size_t g, size_t i, 
                                         const History & history, 
                                         Lattice & L, double T, 
                                         const History & fixed, 
                                         std::vector<std::string> ext) const
{
  return history.subset(ext).derivative<double>().zero();
}

History SlipHardening::d_hist_d_h_ext(const Symmetric & stress,
                                      const Orientation & Q,
                                      const History & history,
                                      Lattice & L, double T, const SlipRule & R,
                                      const History & fixed,
                                      std::vector<std::string> ext) const
{
  return blank_hist().history_derivative(history.subset(ext)).zero();
}

bool SlipHardening::use_nye() const
{
  return false;
}

FixedStrengthHardening::FixedStrengthHardening(ParameterSet & params) :
    SlipHardening(params),
    strengths_(params.get_object_parameter_vector<Interpolate>("strengths"))
{
  init_cache_();
}

std::string FixedStrengthHardening::type()
{
  return "FixedStrengthHardening";
}

std::unique_ptr<NEMLObject> FixedStrengthHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<FixedStrengthHardening>(params);
}

ParameterSet FixedStrengthHardening::parameters()
{
  ParameterSet pset(FixedStrengthHardening::type());

  pset.add_parameter<std::vector<NEMLObject>>("strengths");

  return pset;
}

std::vector<std::string> FixedStrengthHardening::varnames() const
{
  return {};
}

void FixedStrengthHardening::set_varnames(std::vector<std::string> vars)
{
  init_cache_();
}

void FixedStrengthHardening::populate_hist(History & history) const
{
}

void FixedStrengthHardening::init_hist(History & history) const
{
}

double FixedStrengthHardening::hist_to_tau(size_t g, size_t i, 
                                           const History & history,
                                           Lattice & L,
                                           double T, const History & fixed) const
{
  return strengths_[L.flat(g,i)]->value(T);
}

History FixedStrengthHardening::d_hist_to_tau(size_t g, size_t i, 
                                              const History & history,
                                              Lattice & L,
                                              double T, 
                                              const History & fixed) const
{
  return cache(CacheType::DOUBLE);
}

History FixedStrengthHardening::hist(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history, 
                                     Lattice & L, double T, const SlipRule & R, 
                                     const History & fixed) const
{
  return blank_hist();
}

History FixedStrengthHardening::d_hist_d_s(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history,
                                           Lattice & L, double T, 
                                           const SlipRule & R,
                                           const History & fixed) const
{
  return blank_hist().derivative<Symmetric>();
}

History FixedStrengthHardening::d_hist_d_h(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history, 
                                           Lattice & L,
                                           double T, const SlipRule & R, 
                                           const History & fixed) const
{
  return blank_hist().derivative<History>();
}

VocePerSystemHardening::VocePerSystemHardening(ParameterSet & params) :
    SlipHardening(params),
    initial_(params.get_parameter<std::vector<double>>("initial")), 
    k_(params.get_object_parameter_vector<Interpolate>("k")), 
    sat_(params.get_object_parameter_vector<Interpolate>("saturation")), 
    m_(params.get_object_parameter_vector<Interpolate>("m")), 
    varprefix_(params.get_parameter<std::string>("varprefix"))
{
  varnames_.resize(size_());
  for (size_t i = 0; i < size_(); i++) {
    varnames_[i] = varprefix_+std::to_string(i);
  }

  init_cache_();
}

std::string VocePerSystemHardening::type()
{
  return "VocePerSystemHardening";
}

std::unique_ptr<NEMLObject> VocePerSystemHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<VocePerSystemHardening>(params);
}

ParameterSet VocePerSystemHardening::parameters()
{
  ParameterSet pset(VocePerSystemHardening::type());

  pset.add_parameter<std::vector<double>>("initial");
  pset.add_parameter<std::vector<NEMLObject>>("k");
  pset.add_parameter<std::vector<NEMLObject>>("saturation");
  pset.add_parameter<std::vector<NEMLObject>>("m");

  pset.add_optional_parameter<std::string>("varprefix", 
                                           std::string("strength"));

  return pset;
}

std::vector<std::string> VocePerSystemHardening::varnames() const
{
  return varnames_;
}

void VocePerSystemHardening::set_varnames(std::vector<std::string> vars)
{
  if (vars.size() != size_()) 
    throw std::logic_error("New and old varname sizes do not match");
  varnames_ = vars;
  init_cache_();
}

void VocePerSystemHardening::populate_hist(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void VocePerSystemHardening::init_hist(History & history) const
{
  size_t i = 0;
  for (auto vn : varnames_) {
    history.get<double>(vn) = initial_[i];
    i++;
  }
}

double VocePerSystemHardening::hist_to_tau(size_t g, size_t i, 
                                           const History & history,
                                           Lattice & L,
                                           double T, const History & fixed) const
{
  consistency_(L);
  return history.get<double>(varnames_[L.flat(g,i)]);
}

History VocePerSystemHardening::d_hist_to_tau(size_t g, size_t i, 
                                              const History & history,
                                              Lattice & L,
                                              double T, 
                                              const History & fixed) const
{
  History res = cache(CacheType::DOUBLE);
  // This works because the above zeros out the vector
  res.get<double>(varnames_[L.flat(g,i)]) = 1.0;
  return res;
}

History VocePerSystemHardening::hist(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history, 
                                     Lattice & L, double T, const SlipRule & R, 
                                     const History & fixed) const
{
  consistency_(L);
  
  History res = blank_hist();
  
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      size_t k = L.flat(g,i);

      res.get<double>(varnames_[k]) = 
          k_[k]->value(T) * std::pow(
              1.0 - (history.get<double>(varnames_[k]) -initial_[k]) / 
              (sat_[k]->value(T) - initial_[k]), 
              m_[k]->value(T)) * 
          R.slip(g, i, stress, Q, history, L, T, fixed);
    }
  }

  return res;
}

History VocePerSystemHardening::d_hist_d_s(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history,
                                           Lattice & L, double T, 
                                           const SlipRule & R,
                                           const History & fixed) const
{
  History res = blank_hist().derivative<Symmetric>();

  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      size_t k = L.flat(g,i);
      
      res.get<Symmetric>(varnames_[k]) = 
          k_[k]->value(T) * std::pow(
              1.0 - (history.get<double>(varnames_[k]) -initial_[k]) / 
              (sat_[k]->value(T) - initial_[k]), 
              m_[k]->value(T)) * 
          R.d_slip_d_s(g, i, stress, Q, history, L, T, fixed);
    }
  }

  return res;
}

History VocePerSystemHardening::d_hist_d_h(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history, 
                                           Lattice & L,
                                           double T, const SlipRule & R, 
                                           const History & fixed) const
{
  History res = blank_hist().derivative<History>();

  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      size_t k = L.flat(g,i);

      // Self part
      res.get<double>(varnames_[k]+"_"+varnames_[k]) = 
          - k_[k]->value(T) * m_[k]->value(T) / (sat_[k]->value(T) -
                                                 initial_[k]) * 
          std::pow(1.0 - (
                  history.get<double>(varnames_[k]) -initial_[k]) / 
              (sat_[k]->value(T) - initial_[k]), m_[k]->value(T)-1.0) * 
          R.slip(g, i, stress, Q, history, L, T, fixed);

      
      History dslip = R.d_slip_d_h(g, i, stress, Q, history, L, T, fixed);
      for (size_t j = 0; j < L.ntotal(); j++) {
        std::string other = varnames_[j];

        res.get<double>(varnames_[k]+"_"+other) += 
          k_[k]->value(T) * std::pow(
              1.0 - (history.get<double>(varnames_[k]) -initial_[k]) / 
              (sat_[k]->value(T) - initial_[k]), 
              m_[k]->value(T)) * 
            dslip.get<double>(other);
      }
    }
  }

  return res;
}

void VocePerSystemHardening::consistency_(Lattice & L) const
{
  if (L.ntotal() != size_())
    throw std::logic_error("Hardening model size not consistent with lattice!");
}

FASlipHardening::FASlipHardening(ParameterSet & params) :
    SlipHardening(params),
    k_(params.get_object_parameter_vector<Interpolate>("k")),
    sat_(params.get_object_parameter_vector<Interpolate>("saturation")),
    varprefix_(params.get_parameter<std::string>("varprefix"))
{
  varnames_.resize(size_());
  for (size_t i = 0; i < size_(); i++) {
    varnames_[i] = varprefix_+std::to_string(i);
  }

  init_cache_();
}

std::string FASlipHardening::type()
{
  return "FASlipHardening";
}

std::unique_ptr<NEMLObject> FASlipHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<FASlipHardening>(params);
}

ParameterSet FASlipHardening::parameters()
{
  ParameterSet pset(FASlipHardening::type());

  pset.add_parameter<std::vector<NEMLObject>>("k");
  pset.add_parameter<std::vector<NEMLObject>>("saturation");

  pset.add_optional_parameter<std::string>("varprefix", 
                                           std::string("strength"));

  return pset;
}

std::vector<std::string> FASlipHardening::varnames() const
{
  return varnames_;
}

void FASlipHardening::set_varnames(std::vector<std::string> vars)
{
  if (vars.size() != size_()) 
    throw std::logic_error("New and old varname sizes do not match");
  varnames_ = vars;
  init_cache_();
}

void FASlipHardening::populate_hist(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void FASlipHardening::init_hist(History & history) const
{
  size_t i = 0;
  for (auto vn : varnames_) {
    history.get<double>(vn) = 0.0;
    i++;
  }
}

double FASlipHardening::hist_to_tau(size_t g, size_t i, 
                                           const History & history,
                                           Lattice & L,
                                           double T, const History & fixed) const
{
  consistency_(L);
  return history.get<double>(varnames_[L.flat(g,i)]);
}

History FASlipHardening::d_hist_to_tau(size_t g, size_t i, 
                                              const History & history,
                                              Lattice & L,
                                              double T, 
                                              const History & fixed) const
{
  History res = cache(CacheType::DOUBLE);
  // This works because the above zeros out the vector
  res.get<double>(varnames_[L.flat(g,i)]) = 1.0;
  return res;
}

History FASlipHardening::hist(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history, 
                                     Lattice & L, double T, const SlipRule & R, 
                                     const History & fixed) const
{
  consistency_(L);
  
  History res = blank_hist();
  
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      size_t k = L.flat(g,i);
      double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
      std::string vn = varnames_[k];

      res.get<double>(vn) = 
          k_[k]->value(T) * (slip - history.get<double>(vn) / sat_[k]->value(T)
                             * std::fabs(slip));
    }
  }

  return res;
}

History FASlipHardening::d_hist_d_s(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history,
                                           Lattice & L, double T, 
                                           const SlipRule & R,
                                           const History & fixed) const
{
  History res = blank_hist().derivative<Symmetric>();

  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      size_t k = L.flat(g,i);
      double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
      std::string vn = varnames_[k];

      res.get<Symmetric>(vn) = 
          k_[k]->value(T) * (1.0 - history.get<double>(vn) / sat_[k]->value(T)
                             * std::copysign(1.0, slip)) * 
          R.d_slip_d_s(g, i, stress, Q, history, L, T, fixed);
    }
  }

  return res;
}

History FASlipHardening::d_hist_d_h(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history, 
                                           Lattice & L,
                                           double T, const SlipRule & R, 
                                           const History & fixed) const
{
  History res = blank_hist().derivative<History>();

  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      size_t k = L.flat(g,i);
      double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
      std::string vn = varnames_[k];

      // Self part
      res.get<double>(varnames_[k]+"_"+varnames_[k]) = 
          -k_[k]->value(T) / sat_[k]->value(T) * fabs(slip);
      
      History dslip = R.d_slip_d_h(g, i, stress, Q, history, L, T, fixed);
      for (size_t j = 0; j < L.ntotal(); j++) {
        std::string other = varnames_[j];

        res.get<double>(varnames_[k]+"_"+other) += 
            k_[k]->value(T) * (1.0 - history.get<double>(vn) / sat_[k]->value(T)
                               * std::copysign(1.0, slip)) * 
            dslip.get<double>(other);
      }
    }
  }

  return res;
}

void FASlipHardening::consistency_(Lattice & L) const
{
  if (L.ntotal() != size_())
    throw std::logic_error("Hardening model size not consistent with lattice!");
}

GeneralLinearHardening::GeneralLinearHardening(ParameterSet & params) :
    SlipHardening(params),
    M_(params.get_object_parameter<SquareMatrix>("M")), 
    tau_0_(params.get_parameter<std::vector<double>>("tau_0")), 
    absval_(params.get_parameter<bool>("absval")), 
    varprefix_(params.get_parameter<std::string>("varprefix"))
{
  if (M_->n() != tau_0_.size()) {
    throw std::invalid_argument("Hardening matrix and initial strength sizes do not agree!");
  }
  
  varnames_.resize(size());
  for (size_t i = 0; i < size(); i++) {
    varnames_[i] = varprefix_+std::to_string(i);
  }
  init_cache_();
}

std::string GeneralLinearHardening::type()
{
  return "GeneralLinearHardening";
}

std::unique_ptr<NEMLObject> GeneralLinearHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<GeneralLinearHardening>(params);
}

ParameterSet GeneralLinearHardening::parameters()
{
  ParameterSet pset(GeneralLinearHardening::type());

  pset.add_parameter<NEMLObject>("M");
  pset.add_parameter<std::vector<double>>("tau_0");

  pset.add_optional_parameter<bool>("absval", true);
  pset.add_optional_parameter<std::string>("varprefix", 
                                           std::string("strength"));

  return pset;
}

std::vector<std::string> GeneralLinearHardening::varnames() const
{
  return varnames_;
}

void GeneralLinearHardening::set_varnames(std::vector<std::string> vars)
{
  varnames_ = vars;
  init_cache_();
}

void GeneralLinearHardening::populate_hist(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void GeneralLinearHardening::init_hist(History & history) const
{
  size_t i = 0;
  for (auto vn : varnames_) {
    history.get<double>(vn) = tau_0_[i];
    i++;
  }
}

double GeneralLinearHardening::hist_to_tau(size_t g, size_t i, 
                                           const History & history,
                                           Lattice & L,
                                           double T, const History & fixed) const
{
  consistency(L);
  return history.get<double>(varnames_[L.flat(g,i)]) + tau_0_[L.flat(g,i)];
}

History GeneralLinearHardening::d_hist_to_tau(size_t g, size_t i, 
                                              const History & history,
                                              Lattice & L,
                                              double T, 
                                              const History & fixed) const
{
  consistency(L);  
  History res = cache(CacheType::DOUBLE);
  // This works because the above zeros out the vector
  res.get<double>(varnames_[L.flat(g,i)]) = 1.0;
  return res;
}

History GeneralLinearHardening::hist(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history, 
                                     Lattice & L, double T, const SlipRule & R, 
                                     const History & fixed) const
{
  consistency(L); 

  // Vector of slip rates
  FlatVector v(L.ntotal());
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      v.data()[L.flat(g,i)] = R.slip(g, i, stress, Q, history, L, T, fixed);
    }
  }
  if (absval_) {
    for (size_t j = 0; j < L.ntotal(); j++) {
      v.data()[j] = fabs(v.data()[j]);
    }
  }

  // Vector of results
  History res = blank_hist();
  FlatVector resv(L.ntotal(), &(res.get<double>(varnames_[0])));

  // Do the multiplication!
  M_->matvec(v, resv);

  return res;
}

History GeneralLinearHardening::d_hist_d_s(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history,
                                           Lattice & L, double T, 
                                           const SlipRule & R,
                                           const History & fixed) const
{
  consistency(L);
  History res = blank_hist().derivative<Symmetric>();
  
  // Vector of stress derivatives
  std::vector<Symmetric> v(L.ntotal());
  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      v[L.flat(g,i)] = R.d_slip_d_s(g, i, stress, Q, history, L, T, fixed);
    }
  }

  if (absval_) {
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t i = 0; i < L.nslip(g); i++) {
        v[L.flat(g,i)] *= copysign(1.0, R.slip(g, i, stress, Q, history, L, T, fixed));
      }
    }
  }

  // Do the sum
  for (size_t i = 0; i < L.ntotal(); i++) {
    for (size_t j = 0; j < L.ntotal(); j++) {
      res.get<Symmetric>(varnames_[i]) += M_->data()[CINDEX(i,j,L.ntotal())] * v[j];
    }
  }

  return res;
}

History GeneralLinearHardening::d_hist_d_h(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history, 
                                           Lattice & L,
                                           double T, const SlipRule & R, 
                                           const History & fixed) const
{
  consistency(L); 
  auto res = blank_hist().derivative<History>();

  // Do the sum
  for (size_t i = 0; i < L.ntotal(); i++) {
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t l = 0; l < L.nslip(g); l++) {
        size_t j = L.flat(g, l); 
        History curr = R.d_slip_d_h(g, l, stress, Q,  history, L, T, fixed);
        double sign = 1.0;
        if (absval_) {
          sign = copysign(1.0, R.slip(g, l, stress, Q, history, L, T, fixed));
        }
        for (auto vn : varnames_) {
          res.get<double>(varnames_[i]+"_"+vn) += 
              M_->data()[CINDEX(i,j,L.ntotal())] * sign * curr.get<double>(vn);
        }
      }
    }
  }

  return res;
}

History GeneralLinearHardening::d_hist_d_h_ext(const Symmetric & stress, 
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & L, double T, const SlipRule & R,
                                               const History & fixed, 
                                               std::vector<std::string> ext) const
{
  consistency(L);
  History res = blank_hist().history_derivative(history.subset(ext)).zero();

  // Do the sum
  for (size_t i = 0; i < L.ntotal(); i++) {
    for (size_t g = 0; g < L.ngroup(); g++) {
      for (size_t l = 0; l < L.nslip(g); l++) {
        size_t j = L.flat(g, l); 
        History curr = R.d_slip_d_h(g, l, stress, Q,  history, L, T, fixed);
        double sign = 1.0;
        if (absval_) {
          sign = copysign(1.0, R.slip(g, l, stress, Q, history, L, T, fixed));
        }
        for (auto vn : ext) {
          if (curr.contains(vn)) {
            res.get<double>(varnames_[i]+"_"+vn) += 
                M_->data()[CINDEX(i,j,L.ntotal())] * sign * curr.get<double>(vn);
          }
        }
      }
    }
  }

  return res;
}

void GeneralLinearHardening::consistency(Lattice & L) const
{
  if (L.ntotal() != size()) {
    throw std::logic_error("Lattice and hardening matrix sizes do not match");
  }
}

SimpleLinearHardening::SimpleLinearHardening(ParameterSet & params) :
    SlipHardening(params),
    G_(params.get_object_parameter<SquareMatrix>("G")), 
    tau_0_(params.get_parameter<std::vector<double>>("tau_0")), 
    varprefix_(params.get_parameter<std::string>("varprefix"))
{
  if (G_->n() != tau_0_.size()) {
    throw std::invalid_argument("Hardening matrix and initial strength sizes do not agree!");
  }
  
  varnames_.resize(size());
  for (size_t i = 0; i < size(); i++) {
    varnames_[i] = varprefix_+std::to_string(i);
  }
  init_cache_();
}

std::string SimpleLinearHardening::type()
{
  return "SimpleLinearHardening";
}

std::unique_ptr<NEMLObject> SimpleLinearHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<SimpleLinearHardening>(params);
}

ParameterSet SimpleLinearHardening::parameters()
{
  ParameterSet pset(SimpleLinearHardening::type());

  pset.add_parameter<NEMLObject>("G");
  pset.add_parameter<std::vector<double>>("tau_0");

  pset.add_optional_parameter<std::string>("varprefix", 
                                           std::string("slip"));

  return pset;
}

std::vector<std::string> SimpleLinearHardening::varnames() const
{
  return varnames_;
}

void SimpleLinearHardening::set_varnames(std::vector<std::string> vars)
{
  varnames_ = vars;
  init_cache_();
}

void SimpleLinearHardening::populate_hist(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void SimpleLinearHardening::init_hist(History & history) const
{
  size_t i = 0;
  for (auto vn : varnames_) {
    history.get<double>(vn) = 0;
    i++;
  }
}

double SimpleLinearHardening::hist_to_tau(size_t g, size_t i, 
                                           const History & history,
                                           Lattice & L,
                                           double T, const History & fixed) const
{
  consistency(L);

  double v = 0;
  for (size_t j = 0; j < size(); j++) 
    v += (*G_)(i,j) * history.get<double>(varnames_[j]);

  return v + tau_0_[L.flat(g,i)];
}

History SimpleLinearHardening::d_hist_to_tau(size_t g, size_t i, 
                                              const History & history,
                                              Lattice & L,
                                              double T, 
                                              const History & fixed) const
{
  consistency(L);  
  History res = cache(CacheType::DOUBLE);

  for (size_t j = 0; j < size(); j++)
    res.get<double>(varnames_[j]) = (*G_)(i,j);

  return res;
}

History SimpleLinearHardening::hist(const Symmetric & stress, 
                                     const Orientation & Q,
                                     const History & history, 
                                     Lattice & L, double T, const SlipRule & R, 
                                     const History & fixed) const
{
  consistency(L); 

  History res = blank_hist();

  size_t ind = 0;
  for (size_t g = 0; g < L.ngroup(); g++)
    for (size_t i = 0; i < L.nslip(g); i++) {
      res.get<double>(varnames_[ind]) = fabs(R.slip(g, i, stress, Q, history, L,
                                                    T, fixed));
      ind++;
    }

  return res;
}

History SimpleLinearHardening::d_hist_d_s(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history,
                                           Lattice & L, double T, 
                                           const SlipRule & R,
                                           const History & fixed) const
{
  consistency(L);
  History res = blank_hist().derivative<Symmetric>();
 
  size_t ind = 0;
  for (size_t g = 0; g < L.ngroup(); g++)
    for (size_t i = 0; i < L.nslip(g); i++) {
      double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
      res.get<Symmetric>(varnames_[ind]) = copysign(1.0, slip
                                                    ) * R.d_slip_d_s(g, i, stress,
                                                                     Q, history, L,
                                                                     T, fixed);
      ind++;
    }

  return res;
}

History SimpleLinearHardening::d_hist_d_h(const Symmetric & stress, 
                                           const Orientation & Q, 
                                           const History & history, 
                                           Lattice & L,
                                           double T, const SlipRule & R, 
                                           const History & fixed) const
{
  consistency(L); 
  auto res = blank_hist().derivative<History>();

  size_t ind = 0;
  for (size_t g = 0; g < L.ngroup(); g++)
    for (size_t i = 0; i < L.nslip(g); i++) {
      double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
      History dhist = R.d_slip_d_h(g, i, stress, Q, history, L, T, fixed);
      for (size_t j = 0; j < size(); j++)
        res.get<double>(varnames_[ind] + "_" + varnames_[j]) =
            dhist.get<double>(varnames_[j]) * copysign(1.0, slip);
      ind++;
    } 

  return res;
}

History SimpleLinearHardening::d_hist_d_h_ext(const Symmetric & stress, 
                                               const Orientation & Q,
                                               const History & history,
                                               Lattice & L, double T, const SlipRule & R,
                                               const History & fixed, 
                                               std::vector<std::string> ext) const
{
  consistency(L);
  History res = blank_hist().history_derivative(history.subset(ext)).zero();

  size_t ind = 0;
  for (size_t g = 0; g < L.ngroup(); g++)
    for (size_t i = 0; i < L.nslip(g); i++) {
      double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
      History dhist = R.d_slip_d_h(g, i, stress, Q, history, L, T, fixed);
      for (auto vn : ext) {
        if (dhist.contains(vn)) 
          res.get<double>(varnames_[ind] + "_" + vn) = copysign(1.0, slip) *
              dhist.get<double>(vn);
      }
      ind++;
    } 

  return res;
}

void SimpleLinearHardening::consistency(Lattice & L) const
{
  if (L.ntotal() != size()) {
    throw std::logic_error("Lattice and hardening matrix sizes do not match");
  }
}

LANLTiModel::LANLTiModel(ParameterSet & params):
    SlipHardening(params),
    tau_0_(params.get_object_parameter_vector<Interpolate>("tau_0")),
    C_st_(params.get_object_parameter<SquareMatrix>("C_st")),
    mu_(params.get_object_parameter_vector<Interpolate>("mu")),
    k1_(params.get_object_parameter_vector<Interpolate>("k1")), 
    k2_(params.get_object_parameter_vector<Interpolate>("k2")),
    X_s_(params.get_parameter<double>("X_s")),
    inivalue_(params.get_parameter<double>("inivalue")),
    varprefix_(params.get_parameter<std::string>("varprefix")), 
    twinprefix_(params.get_parameter<std::string>("twinprefix"))
{ 

  if (C_st_->n() != nslip_() or C_st_->m() != ntwin_()) {
    throw std::invalid_argument("Twinning interaction matrix and initial strength sizes do not agree!");
  }

  varnames_.resize(size());
  for (size_t i = 0; i < size(); i++) {
    if (i < nslip_()) {
      varnames_[i] = varprefix_+std::to_string(i);
    } 
    else {
      varnames_[i] = twinprefix_+std::to_string(i);	
    }
  }
  init_cache_();
}

std::string LANLTiModel::type()
{
  return "LANLTiModel";
}

std::unique_ptr<NEMLObject> LANLTiModel::initialize(ParameterSet & params)
{
  return neml::make_unique<LANLTiModel>(params);
}

ParameterSet LANLTiModel::parameters()
{
  ParameterSet pset(LANLTiModel::type());

  pset.add_parameter<std::vector<NEMLObject>>("tau_0");
  pset.add_parameter<NEMLObject>("C_st");
  pset.add_parameter<std::vector<NEMLObject>>("mu");
  pset.add_parameter<std::vector<NEMLObject>>("k1");
  pset.add_parameter<std::vector<NEMLObject>>("k2");
  pset.add_optional_parameter<double>("X_s", 0.9);
  pset.add_optional_parameter<double>("inivalue", 0.0);
  pset.add_optional_parameter<std::string>("varprefix", 
                                           std::string("rho"));
  pset.add_optional_parameter<std::string>("twinprefix", 
                                           std::string("slip"));

  return pset;
}

std::vector<std::string> LANLTiModel::varnames() const
{
  return varnames_;
}

void LANLTiModel::set_varnames(std::vector<std::string> vars)
{
  varnames_ = vars;
  init_cache_();
}

void LANLTiModel::populate_hist(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void LANLTiModel::init_hist(History & history) const
{
  for (size_t i = 0; i < size(); i++) {
    if (i < nslip_()) {
      history.get<double>(varnames_[i]) = inivalue_; 
    } 
    else {
      history.get<double>(varnames_[i]) = 0.0;
    }
  }
}

double LANLTiModel::hist_to_tau(size_t g, size_t i, 
                                const History & history,
                                Lattice & L,
                                double T, const History & fixed) const
{
  consistency(L);

  Lattice::SlipType stype = L.slip_type(g,i);	
  if (stype == Lattice::SlipType::Slip) {
    return X_s_ * L.burgers(g,i) * mu_[L.flat(g,i)]->value(T) 
        * history.get<double>(varnames_[L.flat(g,i)])
        + tau_0_[L.flat(g,i)]->value(T);
  } 
  else {
    double v = 0;
    for (size_t g2 = 0; g2 < L.ngroup(); g2++) {
      for (size_t i2 = 0; i2 < L.nslip(g2); i2++) {
        size_t k2 = L.flat(g2,i2);
        Lattice::SlipType otype = L.slip_type(g2,i2);
        if (otype == Lattice::SlipType::Slip) {
          double sqddi = history.get<double>(varnames_[k2]);
          if (sqddi > 0)
            v += (*C_st_)(L.flat(g,i)-nslip_(),k2) * L.burgers(g2,i2)
                * sqddi * sqddi;
        }
      }
    }
    return v * L.burgers(g,i) * mu_[L.flat(g,i)]->value(T)
        + tau_0_[L.flat(g,i)]->value(T);
  }
}  


History LANLTiModel::d_hist_to_tau(size_t g, size_t i, 
                                   const History & history,
                                   Lattice & L,
                                   double T, 
                                   const History & fixed) const
{
  consistency(L);
  History res = cache(CacheType::DOUBLE);

  Lattice::SlipType stype = L.slip_type(g,i);
  if (stype == Lattice::SlipType::Slip) {
    res.get<double>(varnames_[L.flat(g,i)]) = X_s_ * L.burgers(g,i) * mu_[L.flat(g,i)]->value(T);
  }
  else {
    for (size_t g2 = 0; g2 < L.ngroup(); g2++) {
      for (size_t i2 = 0; i2 < L.nslip(g2); i2++) {
        size_t k2 = L.flat(g2,i2);
        Lattice::SlipType otype = L.slip_type(g2,i2);
        if (otype == Lattice::SlipType::Slip) {
          double sqddi = history.get<double>(varnames_[k2]);
          if (sqddi > 0)
            res.get<double>(varnames_[k2]) = 2.0*(*C_st_)(L.flat(g,i)-nslip_(),k2)
                * L.burgers(g2,i2) * mu_[L.flat(g,i)]->value(T) * L.burgers(g,i)
                * sqddi;
        }
      }
    }
  }
  return res;
}

History LANLTiModel::hist(const Symmetric & stress, 
                          const Orientation & Q,
                          const History & history, 
                          Lattice & L, double T, const SlipRule & R, 
                          const History & fixed) const
{
  consistency(L); 

  History res = blank_hist();

  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      size_t k = L.flat(g,i);
      Lattice::SlipType stype = L.slip_type(g,i);
      if (stype == Lattice::SlipType::Slip) {
        res.get<double>(varnames_[k]) = 
            0.5*(k1_[k]->value(T) - k2_[k]->value(T) * history.get<double>(varnames_[k])) * 
            fabs(R.slip(g,i,stress,Q,history,L,T,fixed));
      } 
      else {
        res.get<double>(varnames_[k]) = fabs(R.slip(g, i,
                                                    stress, Q, history, L, T, fixed));  
      }
    }
  }
  return res;
}

History LANLTiModel::d_hist_d_s(const Symmetric & stress, 
                                const Orientation & Q, 
                                const History & history,
                                Lattice & L, double T, 
                                const SlipRule & R,
                                const History & fixed) const
{
  consistency(L);
  History res = blank_hist().derivative<Symmetric>();

  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      Lattice::SlipType stype = L.slip_type(g,i);   
      size_t k = L.flat(g,i);
      if (stype == Lattice::SlipType::Slip) {
        double slip = R.slip(g, i, stress, Q, history, L, T, fixed);        
        res.get<Symmetric>(varnames_[k]) = 
            0.5*(k1_[k]->value(T) - k2_[k]->value(T) * history.get<double>(varnames_[k]))
            * R.d_slip_d_s(g, i, stress, Q, history, L, T, fixed) * copysign(1.0, slip);
      }
      else {
        double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
        res.get<Symmetric>(varnames_[k]) = copysign(1.0, slip)
            * R.d_slip_d_s(g, i, stress, Q, history, L, T, fixed); 								
      }
    }
  }
  return res;
}

History LANLTiModel::d_hist_d_h(const Symmetric & stress, 
                                const Orientation & Q, 
                                const History & history, 
                                Lattice & L,
                                double T, const SlipRule & R, 
                                const History & fixed) const
{
  consistency(L); 
  auto res = blank_hist().derivative<History>();

  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      Lattice::SlipType stype = L.slip_type(g,i);
      size_t k = L.flat(g,i);
      if (stype == Lattice::SlipType::Slip) {	
        History dslip = R.d_slip_d_h(g, i, stress, Q, history, L, T, fixed);
        double slip = R.slip(g, i, stress, Q, history, L, T, fixed);	 
        res.get<double>(varnames_[k] + "_" + varnames_[k]) = -0.5 * k2_[k]->value(T) * std::fabs(slip);
        // other parts		 
        for (size_t j = 0; j < size(); j++) {
          std::string other = varnames_[j];
          res.get<double>(varnames_[k] + "_" + other) += 
              0.5*(k1_[k]->value(T) - k2_[k]->value(T) * history.get<double>(varnames_[k]))
              * dslip.get<double>(other) * copysign(1.0, slip);
        } 
      }
      else {
        History dslip = R.d_slip_d_h(g, i, stress, Q, history, L, T, fixed);
        double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
        for (size_t j = 0; j < size(); j++) {
          std::string other = varnames_[j];
          res.get<double>(varnames_[k] + "_" + other) = dslip.get<double>(other)
              * copysign(1.0, slip);
        }
      }
    }
  }
  return res;
}

History LANLTiModel::d_hist_d_h_ext(const Symmetric & stress, 
                                    const Orientation & Q,
                                    const History & history,
                                    Lattice & L, double T, const SlipRule & R,
                                    const History & fixed, 
                                    std::vector<std::string> ext) const
{
  consistency(L);
  History res = blank_hist().history_derivative(history.subset(ext)).zero();


  for (size_t g = 0; g < L.ngroup(); g++) {
    for (size_t i = 0; i < L.nslip(g); i++) {
      Lattice::SlipType stype = L.slip_type(g,i); 
      size_t k = L.flat(g,i);
      if (stype == Lattice::SlipType::Slip) {
        History dslip = R.d_slip_d_h(g, i, stress, Q, history, L, T, fixed);
        for (auto vn : ext) {
          res.get<double>(varnames_[k] + "_" + vn) = 
              0.5*(k1_[k]->value(T) - k2_[k]->value(T) * history.get<double>(varnames_[k]))
              * dslip.get<double>(vn);
        }				
      }
      else {
        double slip = R.slip(g, i, stress, Q, history, L, T, fixed);
        History dslip = R.d_slip_d_h(g, i, stress, Q, history, L, T, fixed);
        for (auto vn : ext) {
          if (dslip.contains(vn)) { 
            res.get<double>(varnames_[k] + "_" + vn) = 
                dslip.get<double>(vn) * copysign(1.0, slip);
          }
        }
      }
    }
  }
  return res;
}

void LANLTiModel::consistency(Lattice & L) const
{
  if (L.ntotal() != size()) {
    throw std::logic_error("Lattice and hardening matrix sizes do not match");
  }
}

SlipSingleHardening::SlipSingleHardening(ParameterSet & params) :
    SlipHardening(params)
{

}

double SlipSingleHardening::hist_to_tau(
    size_t g, size_t i, const History & history, Lattice & L, double T, 
    const History & fixed) const
{
  return hist_map(history, T, fixed);
}

History SlipSingleHardening::d_hist_to_tau(
    size_t g, size_t i, const History & history, Lattice & L, double T,
    const History & fixed) const
{
  return d_hist_map(history, T, fixed);
}

SlipSingleStrengthHardening::SlipSingleStrengthHardening(ParameterSet & params)
  : 
      SlipSingleHardening(params),
      var_name_(params.get_parameter<std::string>("var_name"))
{

}

std::vector<std::string> SlipSingleStrengthHardening::varnames() const
{
  return {var_name_};
}

void SlipSingleStrengthHardening::set_varnames(std::vector<std::string> vars)
{
  set_variable(vars[0]);
  init_cache_();
}

void SlipSingleStrengthHardening::populate_hist(History & history) const
{
  history.add<double>(var_name_);
}

void SlipSingleStrengthHardening::init_hist(History & history) const
{
  history.get<double>(var_name_) = init_strength();
}

History SlipSingleStrengthHardening::hist(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    Lattice & L, double T, const SlipRule & R, 
    const History & fixed) const
{
  History res = blank_hist();
  res.get<double>(var_name_) = hist_rate(stress, Q, history, L, T, R,
                                         fixed);
  return res;
}

History SlipSingleStrengthHardening::d_hist_d_s(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    Lattice & L, double T,
    const SlipRule & R, const History & fixed) const
{
  History res = blank_hist().derivative<Symmetric>();
  res.get<Symmetric>(var_name_) = d_hist_rate_d_stress(stress,
                                                               Q, 
                                                               history,L, T, R,
                                                               fixed);
  return res;
}

History SlipSingleStrengthHardening::d_hist_d_h(
    const Symmetric & stress, 
    const Orientation & Q,
    const History & history,
    Lattice & L, 
    double T, const SlipRule & R, 
    const History & fixed) const
{
  History res = blank_hist().derivative<History>();
  res.get<double>(var_name_+"_"+var_name_) = d_hist_rate_d_hist(
      stress, Q, history, L, T, R, fixed).get<double>(var_name_);
  return res;
}

History SlipSingleStrengthHardening::d_hist_d_h_ext(
    const Symmetric & stress, const Orientation & Q, const History & history,
    Lattice & L, double T, const SlipRule & R, const History & fixed, 
    std::vector<std::string> ext) const
{
  History res = blank_hist().history_derivative(history.subset(ext)).zero();

  History local_hist = d_hist_rate_d_hist_ext(stress, Q, history,
                                              L, T, R, fixed, ext);

  for (auto name : ext) {
    res.get<double>(var_name_+"_"+name) = local_hist.get<double>(name);
  }

  return res;
}

double SlipSingleStrengthHardening::hist_map(const History & history, 
                                             double T,
                                             const History & fixed) const
{
  return history.get<double>(var_name_) 
      + static_strength(T) 
      + nye_contribution(fixed, T);
}

History SlipSingleStrengthHardening::d_hist_map(const History & history, 
                                                double T, 
                                                const History & fixed) const
{
  History res = cache(CacheType::DOUBLE);
  res.get<double>(var_name_) = 1.0;
  return res;
}

void SlipSingleStrengthHardening::set_variable(std::string name)
{
  var_name_ = name;
}

double SlipSingleStrengthHardening::nye_contribution(const History & fixed,
                                                     double T) const 
{
  if (not use_nye()) return 0;

  if (fixed.contains("nye")) {
    return nye_part(fixed.get<RankTwo>("nye"), T);
  }
  else {
    return 0.0;
  }
}

double SlipSingleStrengthHardening::nye_part(const RankTwo & nye, double T) const
{
  return 0.0;
}

SumSlipSingleStrengthHardening::SumSlipSingleStrengthHardening(ParameterSet &
                                                               params) :
    SlipSingleHardening(params),
    models_(params.get_object_parameter_vector<SlipSingleStrengthHardening>("models"))
{
  for (size_t i = 0; i < nmodels(); i++) {
    models_[i]->set_variable("strength"+std::to_string(i));
  }
  init_cache_();
}

std::vector<std::string> SumSlipSingleStrengthHardening::varnames() const
{
  std::vector<std::string> names;
  
  for (size_t i = 0; i < nmodels(); i++) {
    names.push_back("strength"+std::to_string(i));
  }

  return names;
}

void SumSlipSingleStrengthHardening::set_varnames(std::vector<std::string> vars)
{
  for (size_t i = 0; i < nmodels(); i++) {
    models_[i]->set_variable(vars[i]);
  }
  init_cache_();
}

std::string SumSlipSingleStrengthHardening::type()
{
  return "SumSlipSingleStrengthHardening";
}

std::unique_ptr<NEMLObject> SumSlipSingleStrengthHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<SumSlipSingleStrengthHardening>(params);
}

ParameterSet SumSlipSingleStrengthHardening::parameters()
{
  ParameterSet pset(SumSlipSingleStrengthHardening::type());
  
  pset.add_parameter<std::vector<NEMLObject>>("models");

  return pset;
}

void SumSlipSingleStrengthHardening::populate_hist(History & history) const
{
  for (size_t i = 0; i < nmodels(); i++) {
    history.add<double>("strength"+std::to_string(i));
  }
}

void SumSlipSingleStrengthHardening::init_hist(History & history) const
{
  for (size_t i = 0; i < nmodels(); i++) {
    history.get<double>("strength"+std::to_string(i)) = models_[i]->init_strength();
  }
}

History SumSlipSingleStrengthHardening::hist(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    Lattice & L, double T, const SlipRule & R,
    const History & fixed) const
{
  History res = blank_hist();
  for (size_t i = 0; i < nmodels(); i++) {
    res.get<double>("strength"+std::to_string(i)
                    ) = models_[i]->hist_rate(stress, Q, history, L, T, R, 
                                              fixed);
  }

  return res;
}

History SumSlipSingleStrengthHardening::d_hist_d_s(
    const Symmetric & stress, 
    const Orientation & Q, const History & history,
    Lattice & L, double T,
    const SlipRule & R, const History & fixed) const
{
  History res = blank_hist().derivative<Symmetric>();
  for (size_t i = 0; i < nmodels(); i++) {
    res.get<Symmetric>("strength"+std::to_string(i)) = models_[i]->d_hist_rate_d_stress(stress,
                                                               Q, 
                                                               history,L, T, R,
                                                               fixed);
  }
  return res;
}

History SumSlipSingleStrengthHardening::d_hist_d_h(
    const Symmetric & stress, 
    const Orientation & Q,
    const History & history,
    Lattice & L, 
    double T, const SlipRule & R,
    const History & fixed) const
{
  History res = history.derivative<History>();
  for (size_t i = 0; i < nmodels(); i++) {
    History local_hist = models_[i]->d_hist_rate_d_hist(stress, Q, history, L,
                                                        T, R, fixed);
    for (size_t j = 0; j < nmodels(); j++) {
      res.get<double>("strength"+std::to_string(i)+"_strength"+std::to_string(j)) =
          local_hist.get<double>("strength"+std::to_string(j));
    }
  }
  return res;
}

History SumSlipSingleStrengthHardening::d_hist_d_h_ext(
    const Symmetric & stress, const Orientation & Q, const History & history,
    Lattice & L, double T, const SlipRule & R, const History & fixed, 
    std::vector<std::string> ext) const
{
  History res = blank_hist().history_derivative(history.subset(ext)).zero();

  for (size_t i = 0; i < nmodels(); i++) {
    History local_hist = models_[i]->d_hist_rate_d_hist_ext(stress, Q, history,
                                                            L, T, R, fixed,
                                                            ext);
    for (auto name : ext) {
      res.get<double>("strength"+std::to_string(i)+"_"+name) = 
          local_hist.get<double>(name);
    }
  }

  return res;
}

double SumSlipSingleStrengthHardening::hist_map(const History & history, 
                                                double T,
                                                const History & fixed) const
{
  double sum = 0;
  for (size_t i=0; i < nmodels(); i++) {
    sum += history.get<double>("strength"+std::to_string(i)) +
        models_[i]->static_strength(T) + 
        models_[i]->nye_contribution(fixed, T);
  }
  return sum;
}

History SumSlipSingleStrengthHardening::d_hist_map(const History & history, 
                                                   double T,
                                                   const History & fixed) const
{
  History res = cache(CacheType::DOUBLE);
  for (size_t i=0; i < nmodels(); i++) {
    res.get<double>("strength"+std::to_string(i)) = 1.0;
  }
  return res;
}

size_t SumSlipSingleStrengthHardening::nmodels() const
{
  return models_.size();
}

bool SumSlipSingleStrengthHardening::use_nye() const
{
  for (auto model : models_) {
    if (model->use_nye()) return true;
  }
  return false;
}

PlasticSlipHardening::PlasticSlipHardening(ParameterSet & params) 
  : SlipSingleStrengthHardening(params)
{

}

double PlasticSlipHardening::hist_rate(
    const Symmetric & stress, const Orientation & Q,
    const History & history, Lattice & L, double T, const SlipRule & R,
    const History & fixed) const
{
  double strength = history.get<double>(var_name_);
  
  return hist_factor(strength, L, T, fixed) * R.sum_slip(stress, Q, history, L, T,
                                                  fixed);
}

Symmetric PlasticSlipHardening::d_hist_rate_d_stress(
    const Symmetric & stress, const Orientation & Q, 
    const History & history, Lattice & L, double T,
    const SlipRule & R, const History & fixed) const
{
  double strength = history.get<double>(var_name_);
  return hist_factor(strength, L, T, fixed) * R.d_sum_slip_d_stress(stress, Q,
                                                            history, L, T,
                                                            fixed);
}

History PlasticSlipHardening::d_hist_rate_d_hist(
    const Symmetric & stress, const Orientation & Q, 
    const History & history, Lattice & L, double T,
    const SlipRule & R, const History & fixed) const
{
  double var = history.get<double>(var_name_);
  History dhist = R.d_sum_slip_d_hist(stress, Q, history, L, T, fixed);
  dhist.scalar_multiply(hist_factor(var, L, T, fixed));

  dhist.get<double>(var_name_) += d_hist_factor(var, L, T, fixed
                                                ) * R.sum_slip(
                                                    stress, Q, history, L, T, 
                                                    fixed);

  return dhist;
}

/// Derivative of the scalar law wrt all other scalars
History PlasticSlipHardening::d_hist_rate_d_hist_ext(const Symmetric & stress, 
                                                     const Orientation & Q,
                                                     const History & history,
                                                     Lattice & L, double T,
                                                     const SlipRule & R, 
                                                     const History & fixed, 
                                                     std::vector<std::string> ext) const
{
  History res = history.subset(ext).copy_blank();

  // The only non-zero terms are the derivatives that show up in the
  // d_sum_slip_d_hist
  History dhist = R.d_sum_slip_d_hist(stress, Q, history, L, T, fixed);

  double strength = history.get<double>(var_name_);
  double fact = hist_factor(strength, L, T, fixed);

  for (auto name : ext) {
    if (dhist.contains(name)) {
      res.get<double>(name) = fact * dhist.get<double>(name);
    }
  }

  return res;
}

VoceSlipHardening::VoceSlipHardening(ParameterSet & params) :
    PlasticSlipHardening(params),
    tau_sat_( params.get_object_parameter<Interpolate>("tau_sat")),
    b_(params.get_object_parameter<Interpolate>("b")),
    tau_0_(params.get_object_parameter<Interpolate>("tau_0")), 
    k_(params.get_object_parameter<Interpolate>("k"))
{
  init_cache_();
}

std::string VoceSlipHardening::type()
{
  return "VoceSlipHardening";
}

std::unique_ptr<NEMLObject> VoceSlipHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<VoceSlipHardening>(params);
}

ParameterSet VoceSlipHardening::parameters()
{
  ParameterSet pset(VoceSlipHardening::type());
  
  pset.add_parameter<NEMLObject>("tau_sat");
  pset.add_parameter<NEMLObject>("b");
  pset.add_parameter<NEMLObject>("tau_0");

  pset.add_optional_parameter<NEMLObject>("k", make_constant(0));
  pset.add_optional_parameter<std::string>("var_name", std::string("strength"));

  return pset;
}

double VoceSlipHardening::init_strength() const
{
  return 0.0;
}

double VoceSlipHardening::static_strength(double T) const
{
  return tau_0_->value(T);
}

double VoceSlipHardening::hist_factor(double strength, Lattice & L, 
                                      double T, const History & fixed) const
{
  double tau_sat = tau_sat_->value(T);
  double b = b_->value(T);

  return b * (tau_sat - strength);
}

double VoceSlipHardening::d_hist_factor(double strength, Lattice & L, double T,
                                        const History & fixed) const
{
  double b = b_->value(T);

  return -b;
}

bool VoceSlipHardening::use_nye() const
{
  bool is_constant = true;
  try {
    auto cnst = std::dynamic_pointer_cast<ConstantInterpolate>(k_);
  }
  catch (std::bad_cast & e) {
    is_constant = false;
  }

  if (is_constant and (k_->value(0) == 0.0)) {
    return false;
  }
  return true;
}

double VoceSlipHardening::nye_part(const RankTwo & nye, double T) const
{
  return k_->value(T) * sqrt(nye.norm());
}

LinearSlipHardening::LinearSlipHardening(ParameterSet & params) :
    PlasticSlipHardening(params),
    tau0_(params.get_object_parameter<Interpolate>("tau0")),
    k1_(params.get_object_parameter<Interpolate>("k1")),
    k2_(params.get_object_parameter<Interpolate>("k2"))
{
  init_cache_();
}

std::string LinearSlipHardening::type()
{
  return "LinearSlipHardening";
}

std::unique_ptr<NEMLObject> LinearSlipHardening::initialize(
    ParameterSet & params)
{
  return neml::make_unique<LinearSlipHardening>(params);
}

ParameterSet LinearSlipHardening::parameters()
{
  ParameterSet pset(LinearSlipHardening::type());
  
  pset.add_parameter<NEMLObject>("tau0");
  pset.add_parameter<NEMLObject>("k1");
  pset.add_parameter<NEMLObject>("k2");

  pset.add_optional_parameter<std::string>("var_name", std::string("strength"));

  return pset;
}

double LinearSlipHardening::init_strength() const
{
  return 0.0;
}

double LinearSlipHardening::static_strength(double T) const
{
  return tau0_->value(T);
}

double LinearSlipHardening::hist_factor(double strength, Lattice & L, 
                                      double T, const History & fixed) const
{
  return k1_->value(T);
}

double LinearSlipHardening::d_hist_factor(double strength, Lattice & L, double T,
                                        const History & fixed) const
{
  return 0;
}

bool LinearSlipHardening::use_nye() const
{
  return true;
}

double LinearSlipHardening::nye_part(const RankTwo & nye, double T) const
{
  return k2_->value(T) * nye.norm();
}

} // namespace neml
