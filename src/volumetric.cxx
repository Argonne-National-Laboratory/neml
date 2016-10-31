#include "volumetric.h"
#include "nemlmath.h"

namespace neml {

// Base class implementation
VolumetricModel::VolumetricModel()
{

}

VolumetricModel::~VolumetricModel()
{

}

int VolumetricModel::update(
    const double* const e_inc,
    double T_np1, double T_inc,
    double t_np1, double t_inc,
    double * const h_np1, const double * const h_n,
    double * const s_np1, const double * const s_n,
    double * const A_np1) const
{
  // Get the scalar parameters from the child class implementation
  double e_mean = e_inc[0] + e_inc[1] + e_inc[2];
  double s_mean_n = (s_n[0] + s_n[1] + s_n[2]) / 3.0;
  double s_mean_np1, A_mean;
  int ec = update_mean(e_mean, T_np1, T_inc, t_np1, t_inc, h_np1, h_n, 
                       s_mean_np1, s_mean_n, A_mean);
  
  // ADD! into tensor quantities
  for (size_t i=0; i<3; i++) {
    s_np1[i] += s_mean_np1;
    for (size_t j=0; j<3; j++) {
      A_np1[CINDEX(i,j,6)] += A_mean;
    }
  }

  return ec;
}


// Implementation of ConstantK
VModelConstantK::VModelConstantK(double K) :
    VolumetricModel(), K_(K)
{
  
}

VModelConstantK::~VModelConstantK()
{

}

size_t VModelConstantK::nhist() const
{
  return 0;
}

int VModelConstantK::update_mean(
    double e_inc,
    double T_np1, double T_inc,
    double t_np1, double t_inc,
    double * const h_np1, const double * const h_n,
    double & s_np1, double s_n,
    double & A_np1) const
{
  s_np1 = s_n + K_ * e_inc;
  A_np1 = K_;

  return 0;
}

// Getters
double VModelConstantK::K() const
{
  return K_;
}

}
