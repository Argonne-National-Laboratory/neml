#include "elasticity.h"

#include "nemlmath.h"

#include <algorithm>

namespace neml {

ShearModulus::ShearModulus(double mu) :
    mu_(new ConstantInterpolate(mu))
{

}

ShearModulus::ShearModulus(std::shared_ptr<Interpolate> mu) :
    mu_(mu)
{

}

double ShearModulus::modulus(double T) const
{
  return mu_->value(T);
}



BulkModulus::BulkModulus(std::shared_ptr<Interpolate> K) :
    K_(K)
{
  
}

BulkModulus::BulkModulus(double K) :
    K_(new ConstantInterpolate(K))
{
  
}

double BulkModulus::modulus(double T) const
{
  return K_->value(T);
}

IsotropicLinearElasticModel::IsotropicLinearElasticModel(
      std::shared_ptr<ShearModulus> shear,
      std::shared_ptr<BulkModulus> bulk) :
    shear_(shear), bulk_(bulk)
{

}

int IsotropicLinearElasticModel::C(double T, double * const Cv) const
{
  double G = shear_->modulus(T);
  double K = bulk_->modulus(T);
  double l = K - 2.0/3.0 * G;

  std::fill(Cv, Cv+36, 0.0);

  Cv[0] = 2.0 * G + l;
  Cv[1] = l;
  Cv[2] = l;

  Cv[6] = l;
  Cv[7] = 2.0 * G + l;
  Cv[8] = l;

  Cv[12] = l;
  Cv[13] = l;
  Cv[14] = 2.0 * G + l;

  Cv[21] = 2.0 * G;
  Cv[28] = 2.0 * G;
  Cv[35] = 2.0 * G;

  return 0;
}

int IsotropicLinearElasticModel::S(double T, double * const Sv) const
{
  double G = shear_->modulus(T);
  double K = bulk_->modulus(T);
  double l = K - 2.0/3.0 * G;
  
  double a = (G + l)/(2.0 * G * G + 3 * G * l);
  double b = -l / (4.0 * G * G + 6 * G * l);
  double c = 1.0 / (2.0 * G);

  std::fill(Sv, Sv+36, 0.0);

  Sv[0] = a;
  Sv[1] = b;
  Sv[2] = b;

  Sv[6] = b;
  Sv[7] = a;
  Sv[8] = b;

  Sv[12] = b;
  Sv[13] = b;
  Sv[14] = a;

  Sv[21] = c;
  Sv[28] = c;
  Sv[35] = c;

  return 0;
  
}

} // namespace neml
