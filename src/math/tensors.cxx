#include "tensors.h"

#include "nemlmath.h"

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

Symmetric Symmetric::skew() const
{
  return Symmetric(*this) - trace()/3 * Symmetric::id();
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

SymSym RankFour::to_sym() const
{
  SymSym res;
  
  full2mandel(s_, res.s());

  return res;
}

SymSkew RankFour::to_symskew() const
{
  SymSkew res;

  full2skew(s_, res.s());

  return res;
}

SkewSym RankFour::to_skewsym() const
{
  SkewSym res;

  full2wws(s_, res.s());

  return res;
}

double & SymSym::operator()(size_t i, size_t j)
{
  return s_[i*6+j];
}

const double & SymSym::operator()(size_t i, size_t j) const
{
  return s_[i*6+j];
}

RankFour RankFour::dot(const RankFour & other) const
{
  RankFour res;

  mat_mat(9, 9, 9, this->data(), other.data(), res.s());
  
  return res;
}

RankFour RankFour::dot(const SymSym & other) const
{
  return (*this) * other.to_full();
}

RankFour RankFour::dot(const SymSkew & other) const
{
  return (*this) * other.to_full();
}

RankFour RankFour::dot(const SkewSym & other) const
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

RankFour operator*(const RankFour & a, const SymSym & b)
{
  return a.dot(b);
}

RankFour operator*(const RankFour & a, const SymSkew & b)
{
  return a.dot(b);
}

RankFour operator*(const RankFour & a, const SkewSym & b)
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

/* Start SymSym Tensor */
SymSym::SymSym() :
    Tensor(36)
{
  std::fill(s_, s_+36, 0.0);
}

SymSym::SymSym(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 36) {
    throw std::invalid_argument("Input to SymSym must have size 36!");
  }
}

SymSym::SymSym(const std::vector<std::vector<double>> A) :
    Tensor(36)
{
  if (A.size() != 6) {
    throw std::invalid_argument("SymSym must be initiated with a 6x6 array!");
  }
  for (auto vi : A) {
    if (vi.size() != 6) {
      throw std::invalid_argument("SymSym must be initiated with a 6x6 array!");
    }
  }

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      s_[i*6+j] = A[i][j];
    }
  }
}

SymSym::SymSym(double * v) :
    Tensor(v, 36)
{
}

SymSym::SymSym(const double * v) :
    Tensor(v, 36)
{
}

SymSym SymSym::opposite() const
{
  SymSym cpy(*this);
  cpy.negate_();
  return cpy;
}

SymSym SymSym::operator-() const
{
  return opposite();
}

SymSym & SymSym::operator+=(const SymSym & other)
{
  add_(other);
  return *this;
}

SymSym & SymSym::operator-=(const SymSym & other)
{
  return this->operator+=(-other);
}

RankFour SymSym::to_full() const
{
  RankFour res;

  mandel2full(data(), res.s());

  return res;
}

RankFour SymSym::dot(const RankFour & other) const
{
  return this->to_full() * other;
}

SymSym SymSym::dot(const SymSym & other) const
{
  SymSym res;

  mat_mat(6, 6, 6, this->data(), other.data(), res.s());
  
  return res;
}

RankFour SymSym::dot(const SymSkew & other) const
{
  return this->to_full() * other.to_full();
}

RankFour SymSym::dot(const SkewSym & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SymSym::dot(const RankTwo & other) const
{
  return this->to_full() * other;
}

RankTwo SymSym::dot(const Skew & other) const
{
  return this->to_full() * other;
}

Symmetric SymSym::dot(const Symmetric & other) const
{
  Symmetric res;

  mat_vec(this->data(), 6, other.data(), 6, res.s());

  return res;
}

// Binary operators with scalars
SymSym operator*(double s, const SymSym & v)
{
  SymSym cpy(v);
  cpy *= s;
  return cpy;
}

SymSym operator*(const SymSym & v, double s)
{
  return operator*(s, v);
}

SymSym operator/(const SymSym & v, double s)
{
  SymSym cpy(v);
  cpy /= s;
  return cpy;
}

// Various forms of addition
SymSym operator+(const SymSym & a, const SymSym & b)
{
  SymSym cpy(a);
  cpy += b;
  return cpy;
}

SymSym operator-(const SymSym & a, const SymSym & b)
{
  SymSym cpy(a);
  cpy -= b;
  return cpy;
}

// Various forms of multiplication
RankFour operator*(const SymSym & a, const RankFour & b)
{
  return a.dot(b);
}

SymSym operator*(const SymSym & a, const SymSym & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSym & a, const SymSkew & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSym & a, const SkewSym & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSym & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSym &a, const Skew & b)
{
  return a.dot(b);
}

Symmetric operator*(const SymSym & a, const Symmetric & b)
{
  return a.dot(b);
}

SymSym douter(const Symmetric & a, const Symmetric & b)
{
  SymSym res;
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 6; j++) {
      res(i,j) = a.data()[i] * b.data()[j];
    }
  }

  return res;
}

/// io for SymSym tensors
std::ostream & operator<<(std::ostream & os, const SymSym & v)
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

/* Start SymSkew Tensor */
SymSkew::SymSkew() :
    Tensor(18)
{
}

SymSkew::SymSkew(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 18) {
    throw std::invalid_argument("Input to SymSkew must have size 18!");
  }
}

SymSkew::SymSkew(const std::vector<std::vector<double>> A) :
    Tensor(18)
{
  if (A.size() != 6) {
    throw std::invalid_argument("SymSkew must be initiated with a 6x3 array!");
  }
  for (auto vi : A) {
    if (vi.size() != 3) {
      throw std::invalid_argument("SymSkew must be initiated with a 6x3 array!");
    }
  }

  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < 3; j++) {
      s_[i*3+j] = A[i][j];
    }
  }
}

SymSkew::SymSkew(double * v) :
    Tensor(v, 18)
{
}

SymSkew::SymSkew(const double * v) :
    Tensor(v, 18)
{
}

SymSkew SymSkew::opposite() const
{
  SymSkew cpy(*this);
  cpy.negate_();
  return cpy;
}

SymSkew SymSkew::operator-() const
{
  return opposite();
}

SymSkew & SymSkew::operator+=(const SymSkew & other)
{
  add_(other);
  return *this;
}

SymSkew & SymSkew::operator-=(const SymSkew & other)
{
  return this->operator+=(-other);
}

RankFour SymSkew::to_full() const
{
  RankFour res;

  skew2full(data(), res.s());

  return res;
}

RankFour SymSkew::dot(const RankFour & other) const
{
  return this->to_full() * other;
}

RankFour SymSkew::dot(const SymSym & other) const
{
  RankFour res;

  res = this->to_full() * other.to_full();

  return res;
}

RankFour SymSkew::dot(const SymSkew & other) const
{
  RankFour res;

  res = this->to_full() * other.to_full();

  return res;
}

RankFour SymSkew::dot(const SkewSym & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SymSkew::dot(const RankTwo & other) const
{
  return this->to_full() * other;
}

RankTwo SymSkew::dot(const Skew & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SymSkew::dot(const Symmetric & other) const
{
  RankTwo res;
  
  res = this->to_full() * other.to_full();

  return res;
}

// Binary operators with scalars
SymSkew operator*(double s, const SymSkew & v)
{
  SymSkew cpy(v);
  cpy *= s;
  return cpy;
}

SymSkew operator*(const SymSkew & v, double s)
{
  return operator*(s, v);
}

SymSkew operator/(const SymSkew & v, double s)
{
  SymSkew cpy(v);
  cpy /= s;
  return cpy;
}

// Various forms of addition
SymSkew operator+(const SymSkew & a, const SymSkew & b)
{
  SymSkew cpy(a);
  cpy += b;
  return cpy;
}

SymSkew operator-(const SymSkew & a, const SymSkew & b)
{
  SymSkew cpy(a);
  cpy -= b;
  return cpy;
}

// Various forms of multiplication
RankFour operator*(const SymSkew & a, const RankFour & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSkew & a, const SymSym & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSkew & a, const SymSkew & b)
{
  return a.dot(b);
}

RankFour operator*(const SymSkew & a, const SkewSym & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSkew & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSkew &a, const Skew & b)
{
  return a.dot(b);
}

RankTwo operator*(const SymSkew & a, const Symmetric & b)
{
  return a.dot(b);
}

std::ostream & operator<<(std::ostream & os, const SymSkew & v)
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

/* Start SkewSym Tensor */
SkewSym::SkewSym() :
    Tensor(18)
{
  std::fill(s_, s_+18, 0.0);
}

SkewSym::SkewSym(const std::vector<double> v) :
    Tensor(v)
{
  if (v.size() != 18) {
    throw std::invalid_argument("Input to SkewSym must have size 18!");
  }
}

SkewSym::SkewSym(const std::vector<std::vector<double>> A) :
    Tensor(18)
{
  if (A.size() != 3) {
    throw std::invalid_argument("SkewSym must be initiated with a 3x6 array!");
  }
  for (auto vi : A) {
    if (vi.size() != 6) {
      throw std::invalid_argument("SkewSym must be initiated with a 3x6 array!");
    }
  }

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 6; j++) {
      s_[i*6+j] = A[i][j];
    }
  }
}

SkewSym::SkewSym(double * v) :
    Tensor(v, 18)
{
}

SkewSym::SkewSym(const double * v) :
    Tensor(v, 18)
{
}

SkewSym SkewSym::opposite() const
{
  SkewSym cpy(*this);
  cpy.negate_();
  return cpy;
}

SkewSym SkewSym::operator-() const
{
  return opposite();
}

SkewSym & SkewSym::operator+=(const SkewSym & other)
{
  add_(other);
  return *this;
}

SkewSym & SkewSym::operator-=(const SkewSym & other)
{
  return this->operator+=(-other);
}

RankFour SkewSym::to_full() const
{
  RankFour res;

  wws2full(data(), res.s());

  return res;
}

RankFour SkewSym::dot(const RankFour & other) const
{
  return this->to_full() * other;
}

RankFour SkewSym::dot(const SymSym & other) const
{
  RankFour res;

  res = this->to_full() * other.to_full();

  return res;
}

RankFour SkewSym::dot(const SkewSym & other) const
{
  RankFour res;

  res = this->to_full() * other.to_full();

  return res;
}

RankFour SkewSym::dot(const SymSkew & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SkewSym::dot(const RankTwo & other) const
{
  return this->to_full() * other;
}

RankTwo SkewSym::dot(const Skew & other) const
{
  return this->to_full() * other.to_full();
}

RankTwo SkewSym::dot(const Symmetric & other) const
{
  RankTwo res;
  
  res = this->to_full() * other.to_full();

  return res;
}

// Binary operators with scalars
SkewSym operator*(double s, const SkewSym & v)
{
  SkewSym cpy(v);
  cpy *= s;
  return cpy;
}

SkewSym operator*(const SkewSym & v, double s)
{
  return operator*(s, v);
}

SkewSym operator/(const SkewSym & v, double s)
{
  SkewSym cpy(v);
  cpy /= s;
  return cpy;
}

// Various forms of addition
SkewSym operator+(const SkewSym & a, const SkewSym & b)
{
  SkewSym cpy(a);
  cpy += b;
  return cpy;
}

SkewSym operator-(const SkewSym & a, const SkewSym & b)
{
  SkewSym cpy(a);
  cpy -= b;
  return cpy;
}

// Various forms of multiplication
RankFour operator*(const SkewSym & a, const RankFour & b)
{
  return a.dot(b);
}

RankFour operator*(const SkewSym & a, const SymSym & b)
{
  return a.dot(b);
}

RankFour operator*(const SkewSym & a, const SkewSym & b)
{
  return a.dot(b);
}

RankFour operator*(const SkewSym & a, const SymSkew & b)
{
  return a.dot(b);
}

RankTwo operator*(const SkewSym & a, const RankTwo & b)
{
  return a.dot(b);
}

RankTwo operator*(const SkewSym &a, const Skew & b)
{
  return a.dot(b);
}

RankTwo operator*(const SkewSym & a, const Symmetric & b)
{
  return a.dot(b);
}

SkewSym douter(const Skew & a, const Symmetric & b)
{
  SkewSym res;
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 6; j++) {
      res.s()[i*6+j] = a.data()[i] * b.data()[j];
    }
  }

  return res;  
}

std::ostream & operator<<(std::ostream & os, const SkewSym & v)
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

SymSym SymSymSkew_SkewSymSym(const SymSym & S, const Skew & W)
{
  SymSym res;
  
  SymSymSkewmSkewSymSym(S.data(), W.data(), res.s());

  return res;
}

SymSym SymSkewSym_SkewSymSym(const SkewSym & S, const Symmetric & D)
{
  SymSym res;

  SymSkewSymmSkewSymSym(D.data(), S.data(), res.s());

  return res;
}

SymSkew SpecialSymSymSym(const SymSym & S, const Symmetric & D)
{
  SymSkew res;

  SpecialSymSymSym(D.data(), S.data(), res.s());

  return res;
}

} // namespace neml
