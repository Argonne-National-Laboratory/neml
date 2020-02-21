#include "slipharden.h"

namespace neml {

bool SlipHardening::use_nye() const
{
  return false;
}

double SlipSingleHardening::hist_to_tau(
    size_t g, size_t i, const History & history, double T, 
    const History & fixed) const
{
  return hist_map(history, T, fixed);
}

History SlipSingleHardening::d_hist_to_tau(
    size_t g, size_t i, const History & history, double T,
    const History & fixed) const
{
  return d_hist_map(history, T, fixed);
}

SlipSingleStrengthHardening::SlipSingleStrengthHardening(std::string var_name)
  : var_name_(var_name)
{

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
  History res = history.copy_blank();
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
  History res = history.derivative<Symmetric>();
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
  History res = history.derivative<double>();
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
  History res = history.copy_blank();
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
  History res = history.derivative<Symmetric>();
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
  History res = history.derivative<double>();
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
  catch (std::bad_cast e) {
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
