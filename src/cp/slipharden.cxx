#include "slipharden.h"

#include <stdexcept>

namespace neml {

History SlipHardening::blank_hist() const
{
  History hist;
  populate_history(hist);
  hist.zero();
  return hist;
}

bool SlipHardening::use_nye() const
{
  return false;
}

FixedStrengthHardening::FixedStrengthHardening(
    std::vector<std::shared_ptr<Interpolate>> strengths) :
      strengths_(strengths)
{
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
  return blank_hist().derivative<double>();
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
  History res = blank_hist().derivative<double>();
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
        for (size_t k = 0; k < L.ntotal(); k++) {
          res.get<double>(varnames_[i]+"_"+varnames_[k]) += 
              M_->data()[CINDEX(i,j,L.ntotal())] * sign * curr.get<double>(varnames_[k]);
        }
      }
    }
  }

  return res;
}

size_t GeneralLinearHardening::size() const
{
  return tau_0_.size();
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
  History res = history.derivative<History>();
  res.get<double>(var_name_+"_"+var_name_) = d_hist_rate_d_hist(
      stress, Q, history, L, T, R, fixed).get<double>(var_name_);
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
  History res = blank_hist().derivative<double>();
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
  History res = blank_hist().derivative<double>();
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

VoceSlipHardening::VoceSlipHardening(std::shared_ptr<Interpolate> tau_sat,
                                     std::shared_ptr<Interpolate> b,
                                     std::shared_ptr<Interpolate> tau_0,
                                     std::shared_ptr<Interpolate> k,
                                     std::string var_name) :
    PlasticSlipHardening(var_name), tau_sat_(tau_sat), b_(b), tau_0_(tau_0), 
    k_(k)
{
  
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
