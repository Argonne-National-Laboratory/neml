#include "elasticity.h"

#include "nemlmath.h"

#include <algorithm>

namespace neml {

ConstantShearModulus::ConstantShearModulus(double mu) :
    mu_(mu)
{

}

double ConstantShearModulus::modulus(double T) const
{
  return mu_;
}

double ConstantShearModulus::mu() const
{
  return mu_;
}

PolyShearModulus::PolyShearModulus(const std::vector<double> coefs) :
    coefs_(coefs), n_(coefs.size())
{

}

PolyShearModulus::PolyShearModulus(int n, const double * const coefs) :
    n_(n), coefs_(coefs, coefs + n)
{

}

double PolyShearModulus::modulus(double T) const
{
  return polyval(&coefs_[0], coefs_.size(), T);
}

int PolyShearModulus::n() const
{
  return n_;
}

const std::vector<double> & PolyShearModulus::coefs() const
{
  return coefs_;
}


ConstantBulkModulus::ConstantBulkModulus(double K) :
    K_(K)
{
  
}

double ConstantBulkModulus::modulus(double T) const
{
  return K_;
}

double ConstantBulkModulus::K() const
{
  return K_;
}

PolyBulkModulus::PolyBulkModulus(const std::vector<double> coefs) :
    coefs_(coefs), n_(coefs.size())
{

}

PolyBulkModulus::PolyBulkModulus(int n, const double * const coefs) :
    n_(n), coefs_(coefs, coefs + n)
{

}

double PolyBulkModulus::modulus(double T) const
{
  return polyval(&coefs_[0], coefs_.size(), T);
}

int PolyBulkModulus::n() const
{
  return n_;
}

const std::vector<double> & PolyBulkModulus::coefs() const
{
  return coefs_;
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
