#include "elasticity.h"

#include "math/nemlmath.h"
#include "nemlerror.h"

#include <algorithm>
#include <stdexcept>

namespace neml {

SymSym LinearElasticModel::C(double T) const
{
  SymSym res;
  C(T, res.s());
  return res;
}

SymSym LinearElasticModel::S(double T) const
{
  SymSym res;
  S(T, res.s());
  return res;
}

SymSym LinearElasticModel::C(double T, const Orientation & Q) const
{
  return Q.apply(C(T));
}

SymSym LinearElasticModel::S(double T, const Orientation & Q) const
{
  return Q.apply(S(T));
}

double LinearElasticModel::G(double T) const
{
  return G(T, Orientation::createEulerAngles(0,0,0), 
           Vector({1,0,0}), Vector({0,1,0}));
}

double LinearElasticModel::G(double T, const Orientation & Q, const Vector & b,
                             const Vector & n) const
{
  RankTwo dir = outer(b,n)/(b.norm() * n.norm());
  return dir.contract(C(T,Q).dot(dir));
}

IsotropicLinearElasticModel::IsotropicLinearElasticModel(
      std::shared_ptr<Interpolate> m1,
      std::string m1_type,
      std::shared_ptr<Interpolate> m2,
      std::string m2_type) :
    m1_(m1), m2_(m2), m1_type_(m1_type), m2_type_(m2_type)
{
  if (m1_type_ == m2_type) {
    throw std::invalid_argument("Two distinct elastic constants are required!");
  }
  if (valid_types_.find(m1_type) == valid_types_.end()) {
    throw std::invalid_argument("Unknown elastic constant " + m1_type);
  }
  if (valid_types_.find(m2_type) == valid_types_.end()) {
    throw std::invalid_argument("Unknown elastic constant " + m2_type);
  }
}

std::string IsotropicLinearElasticModel::type()
{
  return "IsotropicLinearElasticModel";
}

ParameterSet IsotropicLinearElasticModel::parameters()
{
  ParameterSet pset(IsotropicLinearElasticModel::type());

  pset.add_parameter<NEMLObject>("m1");
  pset.add_parameter<std::string>("m1_type");
  pset.add_parameter<NEMLObject>("m2");
  pset.add_parameter<std::string>("m2_type");

  return pset;
}

std::unique_ptr<NEMLObject> IsotropicLinearElasticModel::initialize(ParameterSet & params)
{
  return neml::make_unique<IsotropicLinearElasticModel>(
      params.get_object_parameter<Interpolate>("m1"),
      params.get_parameter<std::string>("m1_type"),
      params.get_object_parameter<Interpolate>("m2"),
      params.get_parameter<std::string>("m2_type")
      ); 
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

double IsotropicLinearElasticModel::E(double T) const
{
  double G, K;
  get_GK_(T, G, K);

  return 9.0 * K * G / (3.0 * K + G);
}

double IsotropicLinearElasticModel::nu(double T) const
{
  double G, K;
  get_GK_(T, G, K);

  return (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G));
}

double IsotropicLinearElasticModel::K(double T) const
{
  double G, K;
  get_GK_(T, G, K);

  return K;
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
  double m1 = m1_->value(T);
  double m2 = m2_->value(T);

  if (m1_type_ == "shear" and m2_type_ == "bulk") {
    G = m1;
    K = m2;
  }
  else if (m1_type_ == "bulk" and m2_type_ == "shear") {
    G = m2;
    K = m1;
  }
  else if (m1_type_ == "youngs" and m2_type_ == "poissons") {
    G = m1 / (2.0 * (1.0 + m2));
    K = m1 / (3.0 * (1.0 - 2.0 * m2));
  }
  else if (m1_type_ == "poissons" and m2_type_ == "youngs") {
    G = m2 / (2.0 * (1.0 + m1));
    K = m2 / (3.0 * (1.0 - 2.0 * m1));
  }
  else if (m1_type_ == "youngs" and m2_type_ == "shear") {
    G = m2;
    K = m1 * m2 / (3.0 * (3.0 * m2 - m1));
  }
  else if (m1_type_ == "shear" and m2_type_ == "youngs") {
    G = m1;
    K = m2 * m1 / (3.0 * (3.0 * m1 - m2));
  }
  else if (m1_type_ == "youngs" and m2_type_ == "bulk") {
    G = 3.0 * m2 * m1 / (9.0 * m2 - m1);
    K = m2;
  }
  else if (m1_type_ == "bulk" and m2_type_ == "youngs") {
    G = 3.0 * m1 * m2 / (9.0 * m1 - m2);
    K = m1;
  }
  else if (m1_type_ == "poissons" and m2_type_ == "shear") {
    G = m2;
    K = 2.0 * m2 * (1.0 + m1) / (3.0 * (1.0 - 2.0 * m1));
  }
  else if (m1_type_ == "shear" and m2_type_ == "poissons") {
    G = m1;
    K = 2.0 * m1 * (1.0 + m2) / (3.0 * (1.0 - 2.0 * m2));
  }
  else if (m1_type_ == "poissons" and m2_type_ == "bulk") {
    G = 3.0 * m2 * (1.0 - 2.0 * m1) / (2.0 * (1.0 + m1));
    K = m2;
  }
  else if (m1_type_ == "bulk" and m2_type_ == "poissons") {
    G = 3.0 * m1 * (1.0 - 2.0 * m2) / (2.0 * (1.0 + m2));
    K = m1;
  }
  else {
    throw std::invalid_argument("Unknown combination of elastic properties");
  }
}

} // namespace neml
