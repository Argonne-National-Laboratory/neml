#include "deviatoric.h"

#include "nemlmath.h"

#include <algorithm>

namespace neml {

// Base class implementation
DeviatoricModel::DeviatoricModel()
{

}

DeviatoricModel::~DeviatoricModel()
{

}

// Linear elastic model
LEModel::LEModel(std::shared_ptr<ShearModulus> modulus) :
    modulus_(modulus)
{

}

LEModel::~LEModel()
{

}

size_t LEModel::nhist() const
{
  return 0;
}

int LEModel::init_hist(double* const h) const
{
  return 0;
}

int LEModel::update(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{
  double emean = (e_np1[0] + e_np1[1] + e_np1[2]);
  double edev[6];
  std::copy(e_np1, e_np1+6, edev);
  for (int i=0; i<3; i++) {
    edev[i] -= emean / 3.0;
  }
  std::copy(edev, edev+6, s_np1);
  double mu = modulus_->modulus(T_np1);
  for (int i=0; i<6; i++) {
    s_np1[i] *= 2.0 * mu;
  }

  std::fill(A_np1,A_np1+36,0.0);

  A_np1[CINDEX(0,0,6)] = 4.0 / 3.0 * mu;
  A_np1[CINDEX(0,1,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(0,2,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(1,0,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(1,1,6)] = 4.0 / 3.0 * mu;
  A_np1[CINDEX(1,2,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(2,0,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(2,1,6)] = -2.0 / 3.0 * mu;
  A_np1[CINDEX(2,2,6)] = 4.0 / 3.0 * mu;
  A_np1[CINDEX(3,3,6)] = 2.0 * mu;
  A_np1[CINDEX(4,4,6)] = 2.0 * mu;
  A_np1[CINDEX(5,5,6)] = 2.0 * mu;
}

// Rate independent, associative flow model implementation
RIAFModel::RIAFModel(ShearModulus & modulus, YieldSurface & surface, 
          AssociativeHardening & hardening) :
    modulus_(modulus), surface_(surface), hardening_(hardening)
{

}

RIAFModel::~RIAFModel()
{

}

size_t RIAFModel::nhist() const
{

  return 0;
}

int RIAFModel::init_hist(double* const h) const
{
  return 0;
}


int RIAFModel::update(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{

  return 0;
}


} // namespace neml
