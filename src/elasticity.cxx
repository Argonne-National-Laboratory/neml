#include "elasticity.h"

#include "nemlmath.h"

#include <algorithm>

namespace neml {

Modulus::Modulus(double x) :
    x_(new ConstantInterpolate(x))
{

}

Modulus::Modulus(std::shared_ptr<Interpolate> x) :
    x_(x)
{

}

double Modulus::modulus(double T) const
{
  return x_->value(T);
}


ShearModulus::ShearModulus(double mu) :
    Modulus(mu)
{

}

ShearModulus::ShearModulus(std::shared_ptr<Interpolate> mu) :
    Modulus(mu)
{

}


BulkModulus::BulkModulus(std::shared_ptr<Interpolate> K) :
    Modulus(K)
{
  
}

BulkModulus::BulkModulus(double K) :
    Modulus(K)
{
  
}

YoungsModulus::YoungsModulus(double E) :
    Modulus(E)
{

}

YoungsModulus::YoungsModulus(std::shared_ptr<Interpolate> E) :
    Modulus(E)
{

}


PoissonsRatio::PoissonsRatio(std::shared_ptr<Interpolate> nu) :
    Modulus(nu)
{
  
}

PoissonsRatio::PoissonsRatio(double nu) :
    Modulus(nu)
{
  
}


IsotropicLinearElasticModel::IsotropicLinearElasticModel(
      std::shared_ptr<ShearModulus> shear,
      std::shared_ptr<BulkModulus> bulk) :
    a_(shear), b_(bulk), have_modulii_(GK)
{

}

IsotropicLinearElasticModel::IsotropicLinearElasticModel(
      std::shared_ptr<YoungsModulus> youngs,
      std::shared_ptr<PoissonsRatio> poisson) :
    a_(youngs), b_(poisson), have_modulii_(Ev)
{

}

int IsotropicLinearElasticModel::C(double T, double * const Cv) const
{
  double G, K;
  get_GK_(T, G, K);

  return C_calc_(G, K, Cv);
}

int IsotropicLinearElasticModel::S(double T, double * const Sv) const
{
  double G, K;
  get_GK_(T, G, K);

  return S_calc_(G, K, Sv);
}

int IsotropicLinearElasticModel::C_calc_(double G, double K, double * const Cv) const
{
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

int IsotropicLinearElasticModel::S_calc_(double G, double K, double * const Sv) const
{
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

void IsotropicLinearElasticModel::get_GK_(double T, double & G, double & K) const
{
  double a = a_->modulus(T);
  double b = b_->modulus(T);

  if (have_modulii_ == GK) {
    G = a;
    K = b;
  }
  else if (have_modulii_ == Ev) {
    G = a / (2.0 * (1.0 + b));
    K = a / (3.0 * (1.0 - 2.0 * b));
  }
}

} // namespace neml
