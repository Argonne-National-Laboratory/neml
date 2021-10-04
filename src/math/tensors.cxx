#include "math/tensors.h"

#include "math/nemlmath.h"

#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace neml {

Tensor::Tensor(std::size_t n) :
    n_(n), istore_(true)
{
  s_ = new double [n_];
  std::fill(s_, s_+n_, 0.0);
}

Tensor::Tensor(const Tensor & other) :
    n_(other.n()), istore_(true)
{
  s_ = new double[n_];
  std::copy(other.data(), other.data() + n_, s_);
}

Tensor::Tensor(Tensor && other) :
    n_(other.n()), istore_(other.istore())
{
  if (other.istore()) {
    s_ = new double[n_];
    std::copy(other.data(), other.data() + n_, s_);
  }
  else {
    s_ = other.s();
  }
}

Tensor::Tensor(const std::vector<double> flat) : 
  n_(flat.size()), istore_(true)
{
  s_ = new double [n_];
  std::copy(flat.begin(), flat.end(), s_);
}

Tensor::Tensor(double * flat, size_t n) :
    n_(n), istore_(false)
{
  s_ = flat;
}

Tensor::Tensor(const double * flat, size_t n) :
    n_(n), istore_(false)
{
  s_ = const_cast<double*>(flat);
}

Tensor::~Tensor()
{
  if (istore_) {
    delete [] s_;
  }
  s_ = nullptr;
}

Tensor & Tensor::operator=(const Tensor & rhs) {
  if (n_ != rhs.n()) {
    throw std::invalid_argument(
        "Tensors in assignment operator do not have the same size");
  }

  if (this != &rhs) {
    std::copy(rhs.data(), rhs.data() + rhs.n(), s_);
  }

  return *this;
}

void Tensor::copy_data(const double * const indata)
{
  std::copy(indata, indata + n_, s_);
}

Tensor & Tensor::operator=(Tensor && rhs) {
  if (n_ != rhs.n()) {
    throw std::invalid_argument(
        "Tensors in assignment operator do not have the same size");
  }

  if (rhs.istore()) {
    std::copy(rhs.data(), rhs.data() + n_, s_);
  }
  else {
    s_ = rhs.s();
  }

  return *this;
}

Tensor & Tensor::operator*=(double s) 
{
  for (size_t i=0; i<n_; i++) {
    s_[i] *= s;
  }

  return *this;
}

Tensor & Tensor::operator/=(double s)
{
  return this->operator*=(1.0 / s);
}

void Tensor::add_(const Tensor & other)
{
  for (size_t i = 0; i < n_; i++) {
    s_[i] += other.data()[i];
  }
}

void Tensor::negate_()
{
  for (size_t i = 0; i < n_; i++) {
    s_[i] = -s_[i];
  }
}

bool operator==(const Tensor & a, const Tensor & b)
{
  if (a.n() != b.n()) return false;
  for (size_t i = 0; i < a.n(); i++) {
    if (not isclose(a.data()[i], b.data()[i])) return false;
  }
  return true;
}

bool operator!=(const Tensor & a, const Tensor & b)
{
  return not (a == b);
}

Vector::Vector() :
    Tensor(3)
{
  std::fill(s_, s_+3, 0.0);
}

Vector::Vector(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 3) {
    throw std::invalid_argument("Input to vector must have size 3!");
  }
}

Vector::Vector(double * v) :
    Tensor(v, 3)
{
}

Vector::Vector(const double * v) :
    Tensor(v,3)
{
}

Vector Vector::opposite() const
{
  Vector cpy(*this);
  cpy.negate_();
  return cpy;
}

Vector Vector::operator-() const
{
  return opposite();
}

Vector & Vector::operator+=(const Vector & other)
{
  add_(other);
  return *this;
}

Vector & Vector::operator-=(const Vector & other)
{
  return this->operator+=(-other);
}

double & Vector::operator()(size_t i)
{
  return s_[i];
}

const double & Vector::operator()(size_t i) const
{
  return s_[i];
}

double Vector::dot(const Vector & rhs) const
{
  double s = 0.0;
  for (std::size_t i = 0; i < 3; i++) {
    s += s_[i] * rhs.data()[i];
  }
  return s;
}

RankTwo Vector::outer(const Vector & o) const
{
  RankTwo res;

  for (size_t i=0; i<3; i++) {
    for (size_t j=0; j<3; j++) {
      res(i,j) = (*this)(i) * o(j);
    }
  }
  
  return res;
}

double Vector::norm() const
{
  return sqrt(dot(*this));
}

Vector Vector::cross(const Vector & other) const
{
  const double * a1 = data();
  const double * a2 = other.data();

  Vector t;
  t.s()[0] = a1[1] * a2[2] - a1[2] * a2[1];
  t.s()[1] = a1[2] * a2[0] - a1[0] * a2[2];
  t.s()[2] = a1[0] * a2[1] - a1[1] * a2[0];

  return t;
}

Vector & Vector::normalize()
{
  double n = norm();
  for (int i=0; i<3; i++) {
    s_[i] /= n;
  }
  return *this;
}

Vector operator*(double s, const Vector & v)
{
  Vector cpy(v);
  cpy *= s;
  return cpy;
}

Vector operator*(const Vector & v, double s)
{
  return operator*(s, v);
}

Vector operator/(const Vector & v, double s)
{
  Vector cpy(v);
  cpy /= s;
  return cpy;
}

Vector operator+(const Vector & a, const Vector & b)
{
  Vector cpy(a);
  cpy += b;
  return cpy;
}

Vector operator-(const Vector & a, const Vector & b)
{
  Vector cpy(a);
  cpy -= b;
  return cpy;
}

std::ostream & operator<<(std::ostream & os, const Vector & v)
{
  const double * const d = v.data();
  os << "[" << d[0] << " " << d[1] << " " << d[2] << "]";

  return os;
}

RankTwo outer(const Vector & a, const Vector & b)
{
  return a.outer(b);
}

RankTwo::RankTwo() :
    Tensor(9)
{
  std::fill(s_, s_+9, 0.0);
}

RankTwo::RankTwo(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 9) {
    throw std::invalid_argument("Input to RankTwo must have size 9!");
  }
}

RankTwo::RankTwo(double * v) :
    Tensor(v, 9)
{
}

RankTwo::RankTwo(const double * v) :
    Tensor(v, 9)
{
}

RankTwo::RankTwo(const std::vector<std::vector<double>> A) :
    Tensor(9)
{
  if (A.size() != 3) {
    throw std::invalid_argument("RankTwo must be initiated with a 3x3 array");
  }
  for (auto vi : A) {
    if (vi.size() != 3) {
      throw std::invalid_argument("RankTwo must be initiated with a 3x3 array");
    }
  }

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      s_[i*3+j] = A[i][j];
    }
  }
}

RankTwo::RankTwo(const Symmetric & other) :
    RankTwo(other.to_full())
{

}

RankTwo::RankTwo(const Skew & other) :
    RankTwo(other.to_full())
{

}

RankTwo RankTwo::opposite() const
{
  RankTwo cpy(*this);
  cpy.negate_();
  return cpy;
}

RankTwo RankTwo::operator-() const
{
  return opposite();
}

RankTwo & RankTwo::operator+=(const RankTwo & other)
{
  add_(other);
  return *this;
}

RankTwo & RankTwo::operator-=(const RankTwo & other)
{
  return this->operator+=(-other);
}

RankTwo & RankTwo::operator+=(const Symmetric & other)
{
  // TODO: more efficient
  *this += other.to_full();
  return *this;
}

RankTwo & RankTwo::operator-=(const Symmetric & other)
{
  return this->operator+=(-other);
}

RankTwo & RankTwo::operator+=(const Skew & other)
{
  // TODO: more efficient
  *this += other.to_full();
  return *this;
}

RankTwo & RankTwo::operator-=(const Skew & other)
{
  return this->operator+=(-other);
}

double & RankTwo::operator()(size_t i, size_t j)
{
  return s_[i*3+j];
}

const double & RankTwo::operator()(size_t i, size_t j) const
{
  return s_[i*3+j];
}

RankTwo RankTwo::dot(const RankTwo & other) const
{
  RankTwo res;

  mat_mat(3, 3, 3, this->data(), other.data(), res.s());

  return res;
}

RankTwo RankTwo::dot(const Symmetric & other) const
{
  // TODO: more efficient
  return (*this).dot(other.to_full());
}

RankTwo RankTwo::dot(const Skew & other) const
{
  // TODO: more efficient
  return (*this).dot(other.to_full());
}

Vector RankTwo::dot(const Vector & other) const
{
  Vector res;

  mat_vec(this->data(), 3, other.data(), 3, res.s()); 

  return res;
}

RankTwo RankTwo::inverse() const
{
  RankTwo res(*this);
  invert_mat(res.s(), 3);

  return res;
}

RankTwo RankTwo::transpose() const
{
  RankTwo res;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      res(j,i) = (*this)(i,j);
    }
  }

  return res;
}

double RankTwo::norm() const
{
  return sqrt(this->contract(*this));
}

double RankTwo::contract(const RankTwo & other) const
{
  double sum = 0.0;
  for (int i=0; i<9; i++) sum += s_[i] * other.data()[i];
  return sum;
}

double RankTwo::contract(const Skew & other) const
{
  return this->contract(other.to_full());
}

double RankTwo::contract(const Symmetric & other) const
{
  return this->contract(other.to_full());
}

RankTwo operator*(double s, const RankTwo & v)
{
  RankTwo cpy(v);
  cpy *= s;
  return cpy;
}

RankTwo operator*(const RankTwo & v, double s)
{
  return operator*(s, v);
}

RankTwo operator/(const RankTwo & v, double s)
{
  RankTwo cpy(v);
  cpy /= s;
  return cpy;
}

RankTwo operator+(const RankTwo & a, const RankTwo & b)
{
  RankTwo cpy(a);
  cpy += b;
  return cpy;
}

RankTwo operator-(const RankTwo & a, const RankTwo & b)
{
  RankTwo cpy(a);
  cpy -= b;
  return cpy;
}

RankTwo operator+(const RankTwo & a, const Symmetric & b)
{
  RankTwo cpy(a);
  cpy += b;
  return cpy;
}

RankTwo operator-(const RankTwo & a, const Symmetric & b)
{
  RankTwo cpy(a);
  cpy -= b;
  return cpy;
}

RankTwo operator+(const Symmetric & a, const RankTwo & b)
{
  RankTwo cpy(b);
  cpy += a;
  return cpy;
}

RankTwo operator-(const Symmetric & a, const RankTwo & b)
{
  RankTwo cpy(a);
  cpy -= b;
  return cpy;
}

RankTwo operator+(const RankTwo & a, const Skew & b)
{
  RankTwo cpy(a);
  cpy += b;
  return cpy;
}

RankTwo operator-(const RankTwo & a, const Skew & b)
{
  RankTwo cpy(a);
  cpy -= b;
  return cpy;
}

RankTwo operator+(const Skew & a, const RankTwo & b)
{
  RankTwo cpy(b);
  cpy += a;
  return cpy;
}

RankTwo operator-(const Skew & a, const RankTwo & b)
{
  RankTwo cpy(a);
  cpy -= b;
  return cpy;
}

std::ostream & operator<<(std::ostream & os, const RankTwo & v)
{
  const double * const d = v.data();
  os << "[[" << d[0] << " " << d[1] << " " << d[2] << "]" << std::endl;
  os << " [" << d[3] << " " << d[4] << " " << d[5] << "]" << std::endl;
  os << " [" << d[6] << " " << d[7] << " " << d[8] << "]]" << std::endl;

  return os;
}

Vector operator*(const RankTwo & a, const Vector & b)
{
  return a.dot(b);
}

Vector operator*(const Vector & a, const RankTwo & b)
{
  return b.transpose().dot(a);
}

RankTwo operator*(const RankTwo & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const RankTwo & a, const Symmetric & b)
{
  return a.dot(b);
}

RankTwo operator*(const Symmetric & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const RankTwo & a, const Skew & b)
{
  return a.dot(b);
}

RankTwo operator*(const Skew & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const Skew & a, const Symmetric & b)
{
  return a.dot(b);
}

RankTwo operator*(const Symmetric & a, const Skew & b)
{
  return a.dot(b);
}

Symmetric::Symmetric() :
    Tensor(6)
{
  std::fill(s_, s_+6, 0.0);
}

Symmetric::Symmetric(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 6) {
    throw std::invalid_argument("Input to Symmetric must have size 6!");
  }
}

Symmetric::Symmetric(double * v) :
    Tensor(v, 6)
{
}

Symmetric::Symmetric(const double * v) :
    Tensor(v, 6)
{

}

Symmetric::Symmetric(const RankTwo & other) : 
    Tensor(6)
{
  RankTwo sym = 0.5 * (other + other.transpose());
  s_[0] = sym(0,0);
  s_[1] = sym(1,1);
  s_[2] = sym(2,2);
  s_[3] = sym(1,2) * sqrt(2.0);
  s_[4] = sym(0,2) * sqrt(2.0);
  s_[5] = sym(0,1) * sqrt(2.0);
}

RankTwo Symmetric::to_full() const
{
  RankTwo res;

  res(0,0) = s_[0];
  res(1,1) = s_[1];
  res(2,2) = s_[2];

  res(1,2) = s_[3] / sqrt(2.0);
  res(2,1) = s_[3] / sqrt(2.0);

  res(0,2) = s_[4] / sqrt(2.0);
  res(2,0) = s_[4] / sqrt(2.0);

  res(0,1) = s_[5] / sqrt(2.0);
  res(1,0) = s_[5] / sqrt(2.0);

  return res;
}

Symmetric Symmetric::opposite() const
{
  Symmetric cpy(*this);
  cpy.negate_();
  return cpy;
}

Symmetric Symmetric::operator-() const
{
  return opposite();
}

Symmetric & Symmetric::operator+=(const Symmetric & other)
{
  add_(other);
  return *this;
}

Symmetric & Symmetric::operator-=(const Symmetric & other)
{
  return this->operator+=(-other);
}

Symmetric Symmetric::inverse() const
{
  // TODO: more efficient
  return Symmetric((this->to_full()).inverse());
}

Symmetric Symmetric::transpose() const
{
  return Symmetric(*this);
}

double Symmetric::trace() const
{
  return s_[0] + s_[1] + s_[2];
}

Symmetric Symmetric::dev() const
{
  return Symmetric(*this) - trace()/3 * Symmetric::id();
}

double Symmetric::norm() const
{
  // TODO: more efficient
  return this->to_full().norm();
}

Vector Symmetric::dot(const Vector & other) const
{
  // TODO: more efficient
  return (*this).to_full().dot(other);
}

Symmetric Symmetric::dot(const Symmetric & other) const
{
  // TODO: more efficient
  return Symmetric((*this).to_full().dot(other));
}

RankTwo Symmetric::dot(const RankTwo & other) const
{
  // TODO: more efficient
  return (*this).to_full().dot(other);
}

RankTwo Symmetric::dot(const Skew & other) const
{
  // TODO: more efficient
  return (*this).to_full().dot(other);
}

double Symmetric::contract(const RankTwo & other) const
{
  return other.contract(this->to_full());
}

double Symmetric::contract(const Skew & other) const
{
  return other.contract(this->to_full());
}

double Symmetric::contract(const Symmetric & other) const
{
  double sum = 0.0;
  for (size_t i = 0; i<6; i++) sum += s_[i] * other.data()[i];
  return sum;
}

double & Symmetric::operator()(size_t i)
{
  return s_[i];
}

const double & Symmetric::operator()(size_t i) const
{
  return s_[i];
}

Symmetric operator*(double s, const Symmetric & v)
{
  Symmetric cpy(v);
  cpy *= s;
  return cpy;
}

Symmetric operator*(const Symmetric & v, double s)
{
  return operator*(s, v);
}

Symmetric operator/(const Symmetric & v, double s)
{
  Symmetric cpy(v);
  cpy /= s;
  return cpy;
}

Symmetric operator+(const Symmetric & a, const Symmetric & b)
{
  Symmetric cpy(a);
  cpy += b;
  return cpy;
}

Symmetric operator-(const Symmetric & a, const Symmetric & b)
{
  Symmetric cpy(a);
  cpy -= b;
  return cpy;
}

Vector operator*(const Symmetric & a, const Vector & b)
{
  return a.dot(b);
}

Vector operator*(const Vector & a, const Symmetric & b)
{
  return b.transpose().dot(a);
}

Symmetric operator*(const Symmetric & a, const Symmetric & b)
{
  return a.dot(b);
}

std::ostream & operator<<(std::ostream & os, const Symmetric & v)
{
  os << v.to_full();

  return os;
}

Skew::Skew() :
    Tensor(3)
{
  std::fill(s_, s_+3, 0.0);
}

Skew::Skew(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 3) {
    throw std::invalid_argument("Input to Skew must have size 3!");
  }
}

Skew::Skew(double * v) :
    Tensor(v, 3)
{
}

Skew::Skew(const double * v) :
    Tensor(v, 3)
{
}

Skew::Skew(const RankTwo & other) : 
    Tensor(3)
{
  RankTwo skew = 0.5 * (other - other.transpose());
  s_[0] = -skew(1,2);
  s_[1] = skew(0,2);
  s_[2] = -skew(0,1);
}

RankTwo Skew::to_full() const
{
  RankTwo res;

  res(0,0) = 0;
  res(1,1) = 0;
  res(2,2) = 0;

  res(0,1) = -s_[2];
  res(0,2) = s_[1];

  res(1,0) = s_[2];
  res(1,2) = -s_[0];

  res(2,0) = -s_[1];
  res(2,1) = s_[0];

  return res;
}

Skew Skew::opposite() const
{
  Skew cpy(*this);
  cpy.negate_();
  return cpy;
}

Skew Skew::operator-() const
{
  return opposite();
}

Skew & Skew::operator+=(const Skew & other)
{
  add_(other);
  return *this;
}

Skew & Skew::operator-=(const Skew & other)
{
  return this->operator+=(-other);
}

Skew Skew::transpose() const
{
  return -Skew(*this);
}

Vector Skew::dot(const Vector & other) const
{
  // TODO: more efficient
  return (*this).to_full().dot(other);
}

Skew Skew::dot(const Skew & other) const
{
  // TODO: more efficient
  return Skew((*this).to_full().dot(other));
}

RankTwo Skew::dot(const RankTwo & other) const
{
  // TODO: more efficient
  return (*this).to_full().dot(other);
}

RankTwo Skew::dot(const Symmetric & other) const
{
  // TODO: more efficient
  return (*this).to_full().dot(other);
}

double Skew::contract(const RankTwo & other) const
{
  return other.contract(this->to_full());
}

double Skew::contract(const Skew & other) const
{
  return this->to_full().contract(other.to_full());
}

double Skew::contract(const Symmetric & other) const
{
  return this->contract(other.to_full());
}

Skew operator*(double s, const Skew & v)
{
  Skew cpy(v);
  cpy *= s;
  return cpy;
}

Skew operator*(const Skew & v, double s)
{
  return operator*(s, v);
}

Skew operator/(const Skew & v, double s)
{
  Skew cpy(v);
  cpy /= s;
  return cpy;
}

Skew operator+(const Skew & a, const Skew & b)
{
  Skew cpy(a);
  cpy += b;
  return cpy;
}

Skew operator-(const Skew & a, const Skew & b)
{
  Skew cpy(a);
  cpy -= b;
  return cpy;
}

Vector operator*(const Skew & a, const Vector & b)
{
  return a.dot(b);
}

Vector operator*(const Vector & a, const Skew & b)
{
  return b.transpose().dot(a);
}

Skew operator*(const Skew & a, const Skew & b)
{
  return a.dot(b);
}

std::ostream & operator<<(std::ostream & os, const Skew & v)
{
  os << v.to_full();

  return os;
}

/* Start RankFour Tensor */
RankFour::RankFour() :
    Tensor(81)
{
  std::fill(s_, s_+81, 0.0);
}

RankFour::RankFour(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 81) {
    throw std::invalid_argument("Input to RankFour must have size 81!");
  }
}

RankFour::RankFour(const std::vector<std::vector<std::vector<std::vector<double>>>> A) :
    Tensor(81)
{
  if (A.size() != 3) {
    throw std::invalid_argument("RankFour must be initiated with a 3x3x3x3 array!");
  }
  for (auto vi : A) {
    if (vi.size() != 3) {
      throw std::invalid_argument("RankFour must be initiated with a 6x6 array!");
    }
    for (auto vj : vi) {
      if (vj.size() != 3) {
        throw std::invalid_argument("RankFour must be initiated with a 6x6 array!");
      }
      for (auto vk : vj) {
        if (vk.size() != 3) {
          throw std::invalid_argument("RankFour must be initiated with a 6x6 array!");
        }
      }
    }
  }

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        for (size_t l = 0; l < 3; l++) {
          s_[27*i+9*j+3*k+l] = A[i][j][k][l];
        }
      }
    }
  }
}

RankFour::RankFour(double * v) :
    Tensor(v, 81)
{
}

RankFour::RankFour(const double * v) :
    Tensor(v, 81)
{

}

RankFour RankFour::opposite() const
{
  RankFour cpy(*this);
  cpy.negate_();
  return cpy;
}

RankFour RankFour::operator-() const
{
  return opposite();
}

RankFour & RankFour::operator+=(const RankFour & other)
{
  add_(other);
  return *this;
}

RankFour & RankFour::operator-=(const RankFour & other)
{
  return this->operator+=(-other);
}

double & RankFour::operator()(size_t i, size_t j, size_t k, size_t l)
{
  return s_[i*27+j*9+k*3+l];
}

const double & RankFour::operator()(size_t i, size_t j, size_t k, size_t l) const
{
  return s_[i*27+j*9+k*3+l];
}

SymSymR4 RankFour::to_sym() const
{
  SymSymR4 res;
  
  full2mandel(s_, res.s());

  return res;
}

SymSkewR4 RankFour::to_symskew() const
{
  SymSkewR4 res;

  full2skew(s_, res.s());

  return res;
}

SkewSymR4 RankFour::to_skewsym() const
{
  SkewSymR4 res;

  full2wws(s_, res.s());

  return res;
}

double & SymSymR4::operator()(size_t i, size_t j)
{
  return s_[i*6+j];
}

SymSymR4 SymSymR4::inverse() const
{
  SymSymR4 res(*this);
  invert_mat(res.s(), 6);

  return res;
}

SymSymR4 SymSymR4::transpose() const
{
  SymSymR4 res;
  for (size_t i = 0; i<6; i++) {
    for (size_t j = 0; j<6; j++) {
      res(i,j) = (*this)(j,i);
    }
  }
  return res;
}

const double & SymSymR4::operator()(size_t i, size_t j) const
{
  return s_[i*6+j];
}

RankFour RankFour::dot(const RankFour & other) const
{
  RankFour res;

  mat_mat(9, 9, 9, this->data(), other.data(), res.s());
  
  return res;
}

RankFour RankFour::dot(const SymSymR4 & other) const
{
  return (*this) * other.to_full();
}

RankFour RankFour::dot(const SymSkewR4 & other) const
{
  return (*this) * other.to_full();
}

RankFour RankFour::dot(const SkewSymR4 & other) const
{
  return (*this) * other.to_full();
}

RankTwo RankFour::dot(const RankTwo & other) const
{
  RankTwo res;

  mat_vec(this->data(), 9, other.data(), 9, res.s());

  return res;
}

RankTwo RankFour::dot(const Symmetric & other) const
{
  return dot(other.to_full());
}

RankTwo RankFour::dot(const Skew & other) const
{
  return dot(other.to_full());
}

// Binary operators with scalars
RankFour operator*(double s, const RankFour & v)
{
  RankFour cpy(v);
  cpy *= s;
  return cpy;
}

RankFour operator*(const RankFour & v, double s)
{
  return operator*(s, v);
}

RankFour operator/(const RankFour & v, double s)
{
  RankFour cpy(v);
  cpy /= s;
  return cpy;
}

// Various forms of addition
RankFour operator+(const RankFour & a, const RankFour & b)
{
  RankFour cpy(a);
  cpy += b;
  return cpy;
}

RankFour operator-(const RankFour & a, const RankFour & b)
{
  RankFour cpy(a);
  cpy -= b;
  return cpy;
}

// Various forms of multiplication
RankFour operator*(const RankFour & a, const RankFour & b)
{
  return a.dot(b);
}

RankFour operator*(const RankFour & a, const SymSymR4 & b)
{
  return a.dot(b);
}

RankFour operator*(const RankFour & a, const SymSkewR4 & b)
{
  return a.dot(b);
}

RankFour operator*(const RankFour & a, const SkewSymR4 & b)
{
  return a.dot(b);
}

RankTwo operator*(const RankFour & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const RankFour & a, const Symmetric & b)
{
  return a.dot(b);
}

RankTwo operator*(const RankFour & a, const Skew & b)
{
  return a.dot(b);
}

/// io for RankFour tensors
std::ostream & operator<<(std::ostream & os, const RankFour & v)
{
  const double * const d = v.data();
  for (size_t i = 0; i<9; i++) {
    os << "[";
    for (size_t j = 0; j<9; j++) {
      os << d[i*9+j] << " ";
    }
    os << "]" << std::endl;
  }

  return os;
}

/* Start SymSymR4 Tensor */
SymSymR4::SymSymR4() :
    Tensor(36)
{
  std::fill(s_, s_+36, 0.0);
}

SymSymR4::SymSymR4(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 36) {
    throw std::invalid_argument("Input to SymSymR4 must have size 36!");
  }
}

SymSymR4::SymSymR4(const std::vector<std::vector<double>> A) :
    Tensor(36)
{
  if (A.size() != 6) {
    throw std::invalid_argument("SymSymR4 must be initiated with a 6x6 array!");
  }
  for (auto vi : A) {
    if (vi.size() != 6) {
      throw std::invalid_argument("SymSymR4 must be initiated with a 6x6 array!");
    }
  }

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      s_[i*6+j] = A[i][j];
    }
  }
}

SymSymR4::SymSymR4(double * v) :
    Tensor(v, 36)
{
}

SymSymR4::SymSymR4(const double * v) :
    Tensor(v, 36)
{
}

SymSymR4 SymSymR4::opposite() const
{
  SymSymR4 cpy(*this);
  cpy.negate_();
  return cpy;
}

SymSymR4 SymSymR4::operator-() const
{
  return opposite();
}

SymSymR4 & SymSymR4::operator+=(const SymSymR4 & other)
{
  add_(other);
  return *this;
}

SymSymR4 & SymSymR4::operator-=(const SymSymR4 & other)
{
  return this->operator+=(-other);
}

RankFour SymSymR4::to_full() const
{
  RankFour res;

  mandel2full(data(), res.s());

  return res;
}

RankFour SymSymR4::dot(const RankFour & other) const
{
  return this->to_full() * other;
}

SymSymR4 SymSymR4::dot(const SymSymR4 & other) const
{
  SymSymR4 res;

  mat_mat(6, 6, 6, this->data(), other.data(), res.s());
  
  return res;
}

RankFour SymSymR4::dot(const SymSkewR4 & other) const
{
  return this->to_full() * other.to_full();
}

RankFour SymSymR4::dot(const SkewSymR4 & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SymSymR4::dot(const RankTwo & other) const
{
  return this->to_full() * other;
}

RankTwo SymSymR4::dot(const Skew & other) const
{
  return this->to_full() * other;
}

Symmetric SymSymR4::dot(const Symmetric & other) const
{
  Symmetric res;

  mat_vec(this->data(), 6, other.data(), 6, res.s());

  return res;
}

// Binary operators with scalars
SymSymR4 operator*(double s, const SymSymR4 & v)
{
  SymSymR4 cpy(v);
  cpy *= s;
  return cpy;
}

SymSymR4 operator*(const SymSymR4 & v, double s)
{
  return operator*(s, v);
}

SymSymR4 operator/(const SymSymR4 & v, double s)
{
  SymSymR4 cpy(v);
  cpy /= s;
  return cpy;
}

// Various forms of addition
SymSymR4 operator+(const SymSymR4 & a, const SymSymR4 & b)
{
  SymSymR4 cpy(a);
  cpy += b;
  return cpy;
}

SymSymR4 operator-(const SymSymR4 & a, const SymSymR4 & b)
{
  SymSymR4 cpy(a);
  cpy -= b;
  return cpy;
}

// Various forms of multiplication
RankFour operator*(const SymSymR4 & a, const RankFour & b)
{
  return a.dot(b);
}

SymSymR4 operator*(const SymSymR4 & a, const SymSymR4 & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSymR4 & a, const SymSkewR4 & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSymR4 & a, const SkewSymR4 & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSymR4 & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSymR4 &a, const Skew & b)
{
  return a.dot(b);
}

Symmetric operator*(const SymSymR4 & a, const Symmetric & b)
{
  return a.dot(b);
}

SymSymR4 douter(const Symmetric & a, const Symmetric & b)
{
  SymSymR4 res;
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      res(i,j) = a.data()[i] * b.data()[j];
    }
  }

  return res;
}

/// io for SymSymR4 tensors
std::ostream & operator<<(std::ostream & os, const SymSymR4 & v)
{
  const double * const d = v.data();
  for (size_t i = 0; i<6; i++) {
    os << "[";
    for (size_t j = 0; j<6; j++) {
      os << d[i*6+j] << " ";
    }
    os << "]" << std::endl;
  }

  return os;
}

/* Start SymSkewR4 Tensor */
SymSkewR4::SymSkewR4() :
    Tensor(18)
{
}

SymSkewR4::SymSkewR4(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 18) {
    throw std::invalid_argument("Input to SymSkewR4 must have size 18!");
  }
}

SymSkewR4::SymSkewR4(const std::vector<std::vector<double>> A) :
    Tensor(18)
{
  if (A.size() != 6) {
    throw std::invalid_argument("SymSkewR4 must be initiated with a 6x3 array!");
  }
  for (auto vi : A) {
    if (vi.size() != 3) {
      throw std::invalid_argument("SymSkewR4 must be initiated with a 6x3 array!");
    }
  }

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 3; j++) {
      s_[i*3+j] = A[i][j];
    }
  }
}

SymSkewR4::SymSkewR4(double * v) :
    Tensor(v, 18)
{
}

SymSkewR4::SymSkewR4(const double * v) :
    Tensor(v, 18)
{
}

SymSkewR4 SymSkewR4::opposite() const
{
  SymSkewR4 cpy(*this);
  cpy.negate_();
  return cpy;
}

SymSkewR4 SymSkewR4::operator-() const
{
  return opposite();
}

SymSkewR4 & SymSkewR4::operator+=(const SymSkewR4 & other)
{
  add_(other);
  return *this;
}

SymSkewR4 & SymSkewR4::operator-=(const SymSkewR4 & other)
{
  return this->operator+=(-other);
}

RankFour SymSkewR4::to_full() const
{
  RankFour res;

  skew2full(data(), res.s());

  return res;
}

RankFour SymSkewR4::dot(const RankFour & other) const
{
  return this->to_full() * other;
}

RankFour SymSkewR4::dot(const SymSymR4 & other) const
{
  RankFour res;

  res = this->to_full() * other.to_full();

  return res;
}

RankFour SymSkewR4::dot(const SymSkewR4 & other) const
{
  RankFour res;

  res = this->to_full() * other.to_full();

  return res;
}

RankFour SymSkewR4::dot(const SkewSymR4 & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SymSkewR4::dot(const RankTwo & other) const
{
  return this->to_full() * other;
}

RankTwo SymSkewR4::dot(const Skew & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SymSkewR4::dot(const Symmetric & other) const
{
  RankTwo res;
  
  res = this->to_full() * other.to_full();

  return res;
}

// Binary operators with scalars
SymSkewR4 operator*(double s, const SymSkewR4 & v)
{
  SymSkewR4 cpy(v);
  cpy *= s;
  return cpy;
}

SymSkewR4 operator*(const SymSkewR4 & v, double s)
{
  return operator*(s, v);
}

SymSkewR4 operator/(const SymSkewR4 & v, double s)
{
  SymSkewR4 cpy(v);
  cpy /= s;
  return cpy;
}

// Various forms of addition
SymSkewR4 operator+(const SymSkewR4 & a, const SymSkewR4 & b)
{
  SymSkewR4 cpy(a);
  cpy += b;
  return cpy;
}

SymSkewR4 operator-(const SymSkewR4 & a, const SymSkewR4 & b)
{
  SymSkewR4 cpy(a);
  cpy -= b;
  return cpy;
}

// Various forms of multiplication
RankFour operator*(const SymSkewR4 & a, const RankFour & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSkewR4 & a, const SymSymR4 & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSkewR4 & a, const SymSkewR4 & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSkewR4 & a, const SkewSymR4 & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSkewR4 & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSkewR4 &a, const Skew & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSkewR4 & a, const Symmetric & b)
{
  return a.dot(b);
}

std::ostream & operator<<(std::ostream & os, const SymSkewR4 & v)
{
  const double * const d = v.data();
  for (size_t i = 0; i<6; i++) {
    os << "[";
    for (size_t j = 0; j<3; j++) {
      os << d[i*3+j] << " ";
    }
    os << "]" << std::endl;
  }

  return os;
}

/* Start SkewSymR4 Tensor */
SkewSymR4::SkewSymR4() :
    Tensor(18)
{
  std::fill(s_, s_+18, 0.0);
}

SkewSymR4::SkewSymR4(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 18) {
    throw std::invalid_argument("Input to SkewSymR4 must have size 18!");
  }
}

SkewSymR4::SkewSymR4(const std::vector<std::vector<double>> A) :
    Tensor(18)
{
  if (A.size() != 3) {
    throw std::invalid_argument("SkewSymR4 must be initiated with a 3x6 array!");
  }
  for (auto vi : A) {
    if (vi.size() != 6) {
      throw std::invalid_argument("SkewSymR4 must be initiated with a 3x6 array!");
    }
  }

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 6; j++) {
      s_[i*6+j] = A[i][j];
    }
  }
}

SkewSymR4::SkewSymR4(double * v) :
    Tensor(v, 18)
{
}

SkewSymR4::SkewSymR4(const double * v) :
    Tensor(v, 18)
{
}

SkewSymR4 SkewSymR4::opposite() const
{
  SkewSymR4 cpy(*this);
  cpy.negate_();
  return cpy;
}

SkewSymR4 SkewSymR4::operator-() const
{
  return opposite();
}

SkewSymR4 & SkewSymR4::operator+=(const SkewSymR4 & other)
{
  add_(other);
  return *this;
}

SkewSymR4 & SkewSymR4::operator-=(const SkewSymR4 & other)
{
  return this->operator+=(-other);
}

RankFour SkewSymR4::to_full() const
{
  RankFour res;

  wws2full(data(), res.s());

  return res;
}

RankFour SkewSymR4::dot(const RankFour & other) const
{
  return this->to_full() * other;
}

RankFour SkewSymR4::dot(const SymSymR4 & other) const
{
  RankFour res;

  res = this->to_full() * other.to_full();

  return res;
}

RankFour SkewSymR4::dot(const SkewSymR4 & other) const
{
  RankFour res;

  res = this->to_full() * other.to_full();

  return res;
}

RankFour SkewSymR4::dot(const SymSkewR4 & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SkewSymR4::dot(const RankTwo & other) const
{
  return this->to_full() * other;
}

RankTwo SkewSymR4::dot(const Skew & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SkewSymR4::dot(const Symmetric & other) const
{
  RankTwo res;
  
  res = this->to_full() * other.to_full();

  return res;
}

// Binary operators with scalars
SkewSymR4 operator*(double s, const SkewSymR4 & v)
{
  SkewSymR4 cpy(v);
  cpy *= s;
  return cpy;
}

SkewSymR4 operator*(const SkewSymR4 & v, double s)
{
  return operator*(s, v);
}

SkewSymR4 operator/(const SkewSymR4 & v, double s)
{
  SkewSymR4 cpy(v);
  cpy /= s;
  return cpy;
}

// Various forms of addition
SkewSymR4 operator+(const SkewSymR4 & a, const SkewSymR4 & b)
{
  SkewSymR4 cpy(a);
  cpy += b;
  return cpy;
}

SkewSymR4 operator-(const SkewSymR4 & a, const SkewSymR4 & b)
{
  SkewSymR4 cpy(a);
  cpy -= b;
  return cpy;
}

// Various forms of multiplication
RankFour operator*(const SkewSymR4 & a, const RankFour & b)
{
  return a.dot(b);
}

RankFour operator*(const SkewSymR4 & a, const SymSymR4 & b)
{
  return a.dot(b);
}

RankFour operator*(const SkewSymR4 & a, const SkewSymR4 & b)
{
  return a.dot(b);
}

RankFour operator*(const SkewSymR4 & a, const SymSkewR4 & b)
{
  return a.dot(b);
}

RankTwo operator*(const SkewSymR4 & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const SkewSymR4 &a, const Skew & b)
{
  return a.dot(b);
}

RankTwo operator*(const SkewSymR4 & a, const Symmetric & b)
{
  return a.dot(b);
}

SkewSymR4 douter(const Skew & a, const Symmetric & b)
{
  SkewSymR4 res;
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 6; j++) {
      res.s()[i*6+j] = a.data()[i] * b.data()[j];
    }
  }

  return res;  
}

std::ostream & operator<<(std::ostream & os, const SkewSymR4 & v)
{
  const double * const d = v.data();
  for (size_t i = 0; i<3; i++) {
    os << "[";
    for (size_t j = 0; j<6; j++) {
      os << d[i*6+j] << " ";
    }
    os << "]" << std::endl;
  }

  return os;
}

SymSymR4 SymSymR4Skew_SkewSymR4SymR4(const SymSymR4 & S, const Skew & W)
{
  SymSymR4 res;
  
  SymSymR4SkewmSkewSymR4SymR4(S.data(), W.data(), res.s());

  return res;
}

SymSymR4 SymSkewR4Sym_SkewSymR4SymR4(const SkewSymR4 & S, const Symmetric & D)
{
  SymSymR4 res;

  SymSkewR4SymmSkewSymR4SymR4(D.data(), S.data(), res.s());

  return res;
}

SymSkewR4 SpecialSymSymR4Sym(const SymSymR4 & S, const Symmetric & D)
{
  SymSkewR4 res;

  SpecialSymSymR4Sym(D.data(), S.data(), res.s());

  return res;
}

/* Start SymSymSymR6 Tensor */
SymSymSymR6::SymSymSymR6() :
    Tensor(216)
{
  std::fill(s_, s_+216, 0.0);
}

SymSymSymR6::SymSymSymR6(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 216) {
    throw std::invalid_argument("Input to SymSymSymR6 must have size 216!");
  }
}

SymSymSymR6::SymSymSymR6(const std::vector<std::vector<std::vector<double>>> A) :
    Tensor(216)
{
  if (A.size() != 6) {
    throw std::invalid_argument("SymSymSymR6 must be initiated with a 6x6x6 array!");
  }
  for (auto vi : A) {
    if (vi.size() != 6) {
      throw std::invalid_argument("SymSymSymR6 must be initiated with a 6x6x6 array!");
    }
    for (auto vk : vi) {
      if (vk.size() !=6) {
        throw std::invalid_argument("SymSymSymR6 must be initiated with a 6x6x6 array!");
      }
    }
  }

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      for (size_t k = 0; k < 6; k++) {
        s_[i*36+j*6+k] = A[i][j][k];
      }
    }
  }
}

SymSymSymR6::SymSymSymR6(double * v) :
    Tensor(v, 216)
{
}

SymSymSymR6::SymSymSymR6(const double * v) :
    Tensor(v, 216)
{
}

SymSymSymR6 SymSymSymR6::opposite() const
{
  SymSymSymR6 cpy(*this);
  cpy.negate_();
  return cpy;
}

SymSymSymR6 SymSymSymR6::operator-() const
{
  return opposite();
}

SymSymSymR6 & SymSymSymR6::operator+=(const SymSymSymR6 & other)
{
  add_(other);
  return *this;
}

SymSymSymR6 & SymSymSymR6::operator-=(const SymSymSymR6 & other)
{
  return this->operator+=(-other);
}

double & SymSymSymR6::operator()(size_t i, size_t j, size_t k)
{
  return s_[i*36+j*6+k];
}

const double & SymSymSymR6::operator()(size_t i, size_t j, size_t k) const
{
  return s_[i*36+j*6+k];
}

// Binary operators with scalars
SymSymSymR6 operator*(double s, const SymSymSymR6 & v)
{
  SymSymSymR6 cpy(v);
  cpy *= s;
  return cpy;
}

SymSymSymR6 operator*(const SymSymSymR6 & v, double s)
{
  return operator*(s, v);
}

SymSymSymR6 operator/(const SymSymSymR6 & v, double s)
{
  SymSymSymR6 cpy(v);
  cpy /= s;
  return cpy;
}

// Various forms of addition
SymSymSymR6 operator+(const SymSymSymR6 & a, const SymSymSymR6 & b)
{
  SymSymSymR6 cpy(a);
  cpy += b;
  return cpy;
}

SymSymSymR6 operator-(const SymSymSymR6 & a, const SymSymSymR6 & b)
{
  SymSymSymR6 cpy(a);
  cpy -= b;
  return cpy;
}

SymSymR4 SymSymSymR6::dot_i(const Symmetric & other) const
{
  SymSymR4 res;
  
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      for (size_t k = 0; k < 6; k++) {
        res(j,k) += (*this)(i,j,k) * other(i);
      }
    }
  }

  return res;
}

SymSymR4 SymSymSymR6::dot_j(const Symmetric & other) const
{
  SymSymR4 res;
  
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      for (size_t k = 0; k < 6; k++) {
        res(i,k) += (*this)(i,j,k) * other(j);
      }
    }
  }

  return res;
}

SymSymR4 SymSymSymR6::dot_k(const Symmetric & other) const
{
  SymSymR4 res;
  
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      for (size_t k = 0; k < 6; k++) {
        res(i,j) += (*this)(i,j,k) * other(k);
      }
    }
  }

  return res;
}

SymSymSymR6 SymSymSymR6::middle_dot_after(const SymSymR4 & other) const
{
  SymSymSymR6 res;
  
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      for (size_t k = 0; k < 6; k++) {
        for (size_t l = 0; l < 6; l++) {
          res(i,j,k) += (*this)(i,j,l) * other(j,k);
        }
      }
    }
  }

  return res;
}

SymSymSymR6 SymSymSymR6::middle_dot_before(const SymSymR4 & other) const
{
  SymSymSymR6 res;
  
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      for (size_t k = 0; k < 6; k++) {
        for (size_t l = 0; l < 6; l++) {
          res(i,j,k) += other(i,j) * (*this)(j,k,l);
        }
      }
    }
  }

  return res;
}

// Last index outer product
SymSymSymR6 outer_product_k(const SymSymR4 & A, const Symmetric & B)
{
  SymSymSymR6 res;

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      for (size_t k = 0; k < 6; k++) {
        res(i,j,k) += A(i,j) * B(k);
      }
    }
  }

  return res;
}

} // namespace neml
