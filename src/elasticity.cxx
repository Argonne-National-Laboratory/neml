#include "elasticity.h"

#include "math/nemlmath.h"
#include "nemlerror.h"

#include <algorithm>
#include <stdexcept>

namespace neml {

LinearElasticModel::LinearElasticModel(ParameterSet & params) :
    NEMLObject(params)
{

}

SymSymR4 LinearElasticModel::C(double T) const
{
  SymSymR4 res;
  C(T, res.s());
  return res;
}

SymSymR4 LinearElasticModel::S(double T) const
{
  SymSymR4 res;
  S(T, res.s());
  return res;
}

SymSymR4 LinearElasticModel::C(double T, const Orientation & Q) const
{
  return Q.apply(C(T));
}

SymSymR4 LinearElasticModel::S(double T, const Orientation & Q) const
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
      ParameterSet & params) :
    LinearElasticModel(params),
    m1_(params.get_object_parameter<Interpolate>("m1")), 
    m2_(params.get_object_parameter<Interpolate>("m2")),
    m1_type_(params.get_parameter<std::string>("m1_type")),
    m2_type_(params.get_parameter<std::string>("m2_type"))
{
  if (m1_type_ == m2_type_) {
    throw std::invalid_argument("Two distinct elastic constants are required!");
  }
  if (valid_types_.find(m1_type_) == valid_types_.end()) {
    throw std::invalid_argument("Unknown elastic constant " + m1_type_);
  }
  if (valid_types_.find(m2_type_) == valid_types_.end()) {
    throw std::invalid_argument("Unknown elastic constant " + m2_type_);
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
  return neml::make_unique<IsotropicLinearElasticModel>(params); 
}

void IsotropicLinearElasticModel::C(double T, double * const Cv) const
{
  double G, K;
  get_GK_(T, G, K);

  C_calc_(G, K, Cv);
}

void IsotropicLinearElasticModel::S(double T, double * const Sv) const
{
  double G, K;
  get_GK_(T, G, K);

  S_calc_(G, K, Sv);
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

void IsotropicLinearElasticModel::C_calc_(double G, double K, double * const Cv) const
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
}

void IsotropicLinearElasticModel::S_calc_(double G, double K, double * const Sv) const
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

CubicLinearElasticModel::CubicLinearElasticModel(ParameterSet & params) :
    LinearElasticModel(params),
    M1_(params.get_object_parameter<Interpolate>("m1")), 
    M2_(params.get_object_parameter<Interpolate>("m2")), 
    M3_(params.get_object_parameter<Interpolate>("m3")), 
    method_(params.get_parameter<std::string>("method"))
{
  if ((method_ != "moduli") and (method_ != "components")) {
    throw std::invalid_argument("Unknown initialization method " + method_);
  }
}

std::string CubicLinearElasticModel::type()
{
  return "CubicLinearElasticModel";
}

ParameterSet CubicLinearElasticModel::parameters()
{
  ParameterSet pset(CubicLinearElasticModel::type());

  pset.add_parameter<NEMLObject>("m1");
  pset.add_parameter<NEMLObject>("m2");
  pset.add_parameter<NEMLObject>("m3");
  pset.add_parameter<std::string>("method");

  return pset;
}

std::unique_ptr<NEMLObject> CubicLinearElasticModel::initialize(ParameterSet & params)
{
  return neml::make_unique<CubicLinearElasticModel>(params); 
}

void CubicLinearElasticModel::C(double T, double * const Cv) const
{
  double C1, C2, C3;
  get_components_(T, C1, C2, C3);

  std::fill(Cv, Cv+36, 0.0);

  Cv[0] = C1;
  Cv[1] = C2;
  Cv[2] = C2;

  Cv[6] = C2;
  Cv[7] = C1;
  Cv[8] = C2;

  Cv[12] = C2;
  Cv[13] = C2;
  Cv[14] = C1;

  Cv[21] = C3;
  Cv[28] = C3;
  Cv[35] = C3;
}

void CubicLinearElasticModel::S(double T, double * const Sv) const
{
  C(T, Sv);
  invert_mat(Sv, 6);
}

void CubicLinearElasticModel::get_components_(double T, 
                                              double & C1, double & C2,
                                              double & C3) const
{
  if (method_ == "moduli") {
    double E = M1_->value(T);
    double nu = M2_->value(T);
    double mu = M3_->value(T);

    C1 = E/((1+nu)*(1-2*nu)) * (1-nu);
    C2 = E/((1+nu)*(1-2*nu)) * nu;
    C3 = 2.0 * mu;
  }
  else if (method_ == "components") {
    C1 = M1_->value(T);
    C2 = M2_->value(T);
    C3 = M3_->value(T);
  }
  else {
    throw std::invalid_argument("Invalid method in class internal!");
  }
}

TransverseIsotropicLinearElasticModel::TransverseIsotropicLinearElasticModel(
      ParameterSet & params) :
    LinearElasticModel(params),    
    m1_(params.get_object_parameter<Interpolate>("m1")), 
    m2_(params.get_object_parameter<Interpolate>("m2")),
    m3_(params.get_object_parameter<Interpolate>("m3")),
    m4_(params.get_object_parameter<Interpolate>("m4")),
    m5_(params.get_object_parameter<Interpolate>("m5")),
    method_(params.get_parameter<std::string>("method"))
{
  if ((method_ != "components")) {
    throw std::invalid_argument("Unknown initialization method " + method_);
  }
}

std::string TransverseIsotropicLinearElasticModel::type()
{
  return "TransverseIsotropicLinearElasticModel";
}

ParameterSet TransverseIsotropicLinearElasticModel::parameters()
{
  ParameterSet pset(TransverseIsotropicLinearElasticModel::type());

  pset.add_parameter<NEMLObject>("m1");
  pset.add_parameter<NEMLObject>("m2");
  pset.add_parameter<NEMLObject>("m3");
  pset.add_parameter<NEMLObject>("m4");
  pset.add_parameter<NEMLObject>("m5");
  pset.add_parameter<std::string>("method");

  return pset;
}

std::unique_ptr<NEMLObject> TransverseIsotropicLinearElasticModel::initialize(ParameterSet & params)
{
  return neml::make_unique<TransverseIsotropicLinearElasticModel>(params); 
}

void TransverseIsotropicLinearElasticModel::C(double T, double * const Cv) const
{
  double C11, C33, C12, C13, C44;
  get_components_(T, C11, C33, C12, C13, C44);

  std::fill(Cv, Cv+36, 0.0);
  Cv[0] = C11;
  Cv[1] = C12;
  Cv[2] = C13;

  Cv[6] = C12;
  Cv[7] = C11;
  Cv[8] = C13;

  Cv[12] = C13;
  Cv[13] = C13;
  Cv[14] = C33;

  Cv[21] = C44;
  Cv[28] = C44;
  Cv[35] = (C11-C12)/2;
}

void TransverseIsotropicLinearElasticModel::S(double T, double * const Sv) const
{
  C(T, Sv);
  invert_mat(Sv, 6);
}

void TransverseIsotropicLinearElasticModel::get_components_(double T, 
                                              double & C11, double & C33,
                                              double & C12, double & C13,
                                              double & C44) const
{
  if (method_ == "components") {
    C11 = m1_->value(T);
    C33 = m2_->value(T);
    C12 = m3_->value(T);
    C13 = m4_->value(T);
    C44 = m5_->value(T);
  }
  else {
    throw std::invalid_argument("Invalid method in class internal!");
  }
}

} // namespace neml
