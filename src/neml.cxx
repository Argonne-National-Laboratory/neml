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
NEMLModel_ldF::NEMLModel_ldF()
{

}

NEMLModel_ldF::~NEMLModel_ldF()
{

}

size_t NEMLModel_ldF::nstore()
{
  return nhist() + 0;
}

int NEMLModel_ldF::update_ldI(
    const double* const l_inc,
    double T_np1, double T_inc,
    double t_np1, double t_inc,
    double * const h_np1, const double * const h_n,
    double * const s_np1, const double * const s_n,
    double * const A_np1)
{
  assert(false); // Pass on implementing for now
}

int NEMLModel_ldF::update_sd(
     const double* const e_inc,
     double T_np1, double T_inc,
     double t_np1, double t_inc,
     double * const h_np1, const double * const h_n,
     double * const s_np1, const double * const s_n,
     double * const A_np1)
{
  assert(false); // Pass on implementing for now
}

// NEMLModel_ldI implementation
NEMLModel_ldI::NEMLModel_ldI()
{

}

NEMLModel_ldI::~NEMLModel_ldI()
{

}

size_t NEMLModel_ldI::nstore()
{
  return nhist() + 0;
}

int NEMLModel_ldI::update_ldF(
    const double * const F_np1, const double * const F_n,
    double T_np1, double T_inc,
    double t_np1, double t_inc,
    double * const h_np1, const double * const h_n,
    double * const s_np1, const double * const s_n,
    double * const A_np1)
{
  assert(false); // Pass for now
}

int NEMLModel_ldI::update_sd(
    const double* const e_inc,
    double T_np1, double T_inc,
    double t_np1, double t_inc,
    double * const h_np1, const double * const h_n,
    double * const s_np1, const double * const s_n,
    double * const A_np1)
{
  assert(false); // Pass for now
}


// NEMLModel_sd implementation
NEMLModel_sd::NEMLModel_sd()
{

}

NEMLModel_sd::~NEMLModel_sd()
{

}

size_t NEMLModel_sd::nstore()
{
  return nhist() + 4; // The history variables plus a rotation
}

int NEMLModel_sd::update_ldF(
    const double * const F_np1, const double * const F_n,
    double T_np1, double T_inc,
    double t_np1, double t_inc,
    double * const h_np1, const double * const h_n,
    double * const s_np1, const double * const s_n,
    double * const A_np1)
{
  assert(false); // Pass for now
}

int NEMLModel_sd::update_ldI(
    const double* const l_inc,
    double T_np1, double T_inc,
    double t_np1, double t_inc,
    double * const h_np1, const double * const h_n,
    double * const s_np1, const double * const s_n,
    double * const A_np1)
{
  assert(false); // Pass for now
}


// SplitModel_sd implementation
SplitModel_sd::SplitModel_sd()
{

}

SplitModel_sd::~SplitModel_sd()
{

}

size_t SplitModel_sd::nhist()
{

}

int SplitModel_sd::update_sd(
   const double* const e_inc,
   double T_np1, double T_inc,
   double t_np1, double t_inc,
   double * const h_np1, const double * const h_n,
   double * const s_np1, const double * const s_n,
   double * const A_np1)
{

}

} // namespace neml
