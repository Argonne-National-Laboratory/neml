#include "slipharden.h"

#include <stdexcept>

namespace neml {

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
  populate_history(*blank_);
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

FixedStrengthHardening::FixedStrengthHardening(
    std::vector<std::shared_ptr<Interpolate>> strengths) :
      strengths_(strengths)
{
  init_cache_();
}

std::string FixedStrengthHardening::type()
{
  return "FixedStrengthHardening";
}

std::unique_ptr<NEMLObject> FixedStrengthHardening::initialize(ParameterSet & params)
{
  return neml::make_unique<FixedStrengthHardening>(
      params.get_object_parameter_vector<Interpolate>("strengths"));
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

void FixedStrengthHardening::populate_history(History & history) const
{
}

void FixedStrengthHardening::init_history(History & history) const
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

VocePerSystemHardening::VocePerSystemHardening(
    std::vector<double> initial,
    std::vector<std::shared_ptr<Interpolate>> k,
    std::vector<std::shared_ptr<Interpolate>> sat,
    std::vector<std::shared_ptr<Interpolate>> m,
    std::string varprefix) :
      initial_(initial), k_(k), sat_(sat), m_(m), varprefix_(varprefix)
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
  return neml::make_unique<VocePerSystemHardening>(
      params.get_parameter<std::vector<double>>("initial"),
      params.get_object_parameter_vector<Interpolate>("k"),
      params.get_object_parameter_vector<Interpolate>("saturation"),
      params.get_object_parameter_vector<Interpolate>("m"),
      params.get_parameter<std::string>("varprefix"));
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

void VocePerSystemHardening::populate_history(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void VocePerSystemHardening::init_history(History & history) const
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

FASlipHardening::FASlipHardening(
    std::vector<std::shared_ptr<Interpolate>> k,
    std::vector<std::shared_ptr<Interpolate>> sat,
    std::string varprefix) :
      k_(k), sat_(sat), varprefix_(varprefix)
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
  return neml::make_unique<FASlipHardening>(
      params.get_object_parameter_vector<Interpolate>("k"),
      params.get_object_parameter_vector<Interpolate>("saturation"),
      params.get_parameter<std::string>("varprefix"));
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

void FASlipHardening::populate_history(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void FASlipHardening::init_history(History & history) const
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

GeneralLinearHardening::GeneralLinearHardening(std::shared_ptr<SquareMatrix> M, 
                                               std::vector<double> tau_0,
                                               bool absval,
                                               std::string varprefix) :
    M_(M), tau_0_(tau_0), absval_(absval), varprefix_(varprefix)
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
  return neml::make_unique<GeneralLinearHardening>(
      params.get_object_parameter<SquareMatrix>("M"),
      params.get_parameter<std::vector<double>>("tau_0"),
      params.get_parameter<bool>("absval"),
      params.get_parameter<std::string>("varprefix"));
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

void GeneralLinearHardening::populate_history(History & history) const
{
  for (auto vn : varnames_) {
    history.add<double>(vn);
  }
}

void GeneralLinearHardening::init_history(History & history) const
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

SlipSingleStrengthHardening::SlipSingleStrengthHardening(std::string var_name)
  : var_name_(var_name)
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

void SlipSingleStrengthHardening::populate_history(History & history) const
{
  history.add<double>(var_name_);
}

void SlipSingleStrengthHardening::init_history(History & history) const
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

SumSlipSingleStrengthHardening::SumSlipSingleStrengthHardening(
    std::vector<std::shared_ptr<SlipSingleStrengthHardening>> models)
  :   models_(models)
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
  return neml::make_unique<SumSlipSingleStrengthHardening>(
      params.get_object_parameter_vector<SlipSingleStrengthHardening>("models"));
}

ParameterSet SumSlipSingleStrengthHardening::parameters()
{
  ParameterSet pset(SumSlipSingleStrengthHardening::type());
  
  pset.add_parameter<std::vector<NEMLObject>>("models");

  return pset;
}

void SumSlipSingleStrengthHardening::populate_history(History & history) const
{
  for (size_t i = 0; i < nmodels(); i++) {
    history.add<double>("strength"+std::to_string(i));
  }
}

void SumSlipSingleStrengthHardening::init_history(History & history) const
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

PlasticSlipHardening::PlasticSlipHardening(std::string var_name) 
  : SlipSingleStrengthHardening(var_name)
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

VoceSlipHardening::VoceSlipHardening(std::shared_ptr<Interpolate> tau_sat,
                                     std::shared_ptr<Interpolate> b,
                                     std::shared_ptr<Interpolate> tau_0,
                                     std::shared_ptr<Interpolate> k,
                                     std::string var_name) :
    PlasticSlipHardening(var_name), tau_sat_(tau_sat), b_(b), tau_0_(tau_0), 
    k_(k)
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
  return neml::make_unique<VoceSlipHardening>(
      params.get_object_parameter<Interpolate>("tau_sat"),
      params.get_object_parameter<Interpolate>("b"),
      params.get_object_parameter<Interpolate>("tau_0"),
      params.get_object_parameter<Interpolate>("k"));
}

ParameterSet VoceSlipHardening::parameters()
{
  ParameterSet pset(VoceSlipHardening::type());
  
  pset.add_parameter<NEMLObject>("tau_sat");
  pset.add_parameter<NEMLObject>("b");
  pset.add_parameter<NEMLObject>("tau_0");

  pset.add_optional_parameter<NEMLObject>("k", 
                                    std::make_shared<ConstantInterpolate>(0));

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

LinearSlipHardening::LinearSlipHardening(std::shared_ptr<Interpolate> tau0,
                                         std::shared_ptr<Interpolate> k1,
                                         std::shared_ptr<Interpolate> k2,
                                         std::string var_name) :
    PlasticSlipHardening(var_name), tau0_(tau0), k1_(k1), k2_(k2)
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
  return neml::make_unique<LinearSlipHardening>(
      params.get_object_parameter<Interpolate>("tau0"),
      params.get_object_parameter<Interpolate>("k1"),
      params.get_object_parameter<Interpolate>("k2"));
}

ParameterSet LinearSlipHardening::parameters()
{
  ParameterSet pset(LinearSlipHardening::type());
  
  pset.add_parameter<NEMLObject>("tau0");
  pset.add_parameter<NEMLObject>("k1");
  pset.add_parameter<NEMLObject>("k2");

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
