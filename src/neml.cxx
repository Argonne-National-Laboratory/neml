#include "neml.h"

#include <cassert>

namespace neml {

// Base class NEMLModel

NEMLModel::NEMLModel()
{

}

NEMLModel::~NEMLModel()
{

}



// NEMLModel_ldF implementation
NEMLModel_ldF::NEMLModel_ldF() :
    NEMLModel()
{

}

NEMLModel_ldF::~NEMLModel_ldF()
{

}

size_t NEMLModel_ldF::nstore() const
{
  return nhist() + 0;
}

int NEMLModel_ldF::init_store(double * const store) const
{
  assert(false); 
}

int NEMLModel_ldF::update_ldI(
    const double * const l_inc,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{
  assert(false); // Pass on implementing for now
}

int NEMLModel_ldF::update_sd(
     const double * const e_np1, const double * const e_n,
     double T_np1, double T_n,
     double t_np1, double t_n,
     double * const s_np1, const double * const s_n,
     double * const h_np1, const double * const h_n,
     double * const A_np1) const
{
  assert(false); // Pass on implementing for now
}

// NEMLModel_ldI implementation
NEMLModel_ldI::NEMLModel_ldI() :
    NEMLModel()
{

}

NEMLModel_ldI::~NEMLModel_ldI()
{

}

size_t NEMLModel_ldI::nstore() const
{
  return nhist() + 0;
}

int NEMLModel_ldI::init_store(double * const store) const
{
  assert(false); 
}

int NEMLModel_ldI::update_ldF(
    const double * const F_np1, const double * const F_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{
  assert(false); // Pass for now
}

int NEMLModel_ldI::update_sd(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{
  assert(false); // Pass for now
}


// NEMLModel_sd implementation
NEMLModel_sd::NEMLModel_sd() :
    NEMLModel()
{

}

NEMLModel_sd::~NEMLModel_sd()
{

}

size_t NEMLModel_sd::nstore() const
{
  return nhist(); // Get to this later...
}

int NEMLModel_sd::init_store(double * const store) const
{
  init_hist(store); // Get to this later
}

int NEMLModel_sd::update_ldF(
    const double * const F_np1, const double * const F_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{
  assert(false); // Pass for now
}

int NEMLModel_sd::update_ldI(
    const double * const l_inc,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{
  assert(false); // Pass for now
}


// SplitModel_sd implementation
SplitModel_sd::SplitModel_sd(std::shared_ptr<VolumetricModel> vol_model, 
                             std::shared_ptr<DeviatoricModel> dev_model) :
    NEMLModel_sd(), vol_model_(vol_model), dev_model_(dev_model)
{

}

SplitModel_sd::~SplitModel_sd()
{

}

size_t SplitModel_sd::nhist() const
{
  return dev_model_->nhist() + vol_model_->nhist();
}

int SplitModel_sd::init_hist(double * const hist) const
{
  int ier1 = dev_model_->init_hist(hist);
  int ier2 = vol_model_->init_hist(&hist[dev_model_->nhist()]);

  if (ier1 != 0) return ier1;
  return ier2;
}

int SplitModel_sd::update_sd(
   const double * const e_np1, const double * const e_n,
   double T_np1, double T_n,
   double t_np1, double t_n,
   double * const s_np1, const double * const s_n,
   double * const h_np1, const double * const h_n,
   double * const A_np1) const
{
  // In my somewhat confusing convention the volumetric model automatically
  // adds into the existing, nonzero arrays
  dev_model_->update(e_np1, e_n, T_np1, T_n, t_np1, t_n, 
                     s_np1, s_n, h_np1, h_n, A_np1);
  vol_model_->update(e_np1, e_n, T_np1, T_n, t_np1, t_n,
                     s_np1, s_n, 
                     &h_np1[dev_model_->nhist()],
                     &h_n[dev_model_->nhist()],
                     A_np1);
}

} // namespace neml
