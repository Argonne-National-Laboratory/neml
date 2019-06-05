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
}

Tensor::Tensor(const Tensor & other) :
    n_(other.n()), istore_(other.istore())
{
  if (other.istore()) {
    s_ = new double[n_];
    std::copy(other.data(), other.data() + n_, s_);
  }
  else {
    s_ = const_cast<double*>(other.data());
  }
}

Tensor::Tensor(const Tensor && other) :
    n_(other.n()), istore_(other.istore())
{
  if (other.istore()) {
    s_ = new double[n_];
    std::copy(other.data(), other.data() + n_, s_);
  }
  else {
    s_ = const_cast<double*>(other.data());
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

Tensor & Tensor::operator=(Tensor && rhs) {
  if (n_ != rhs.n()) {
    throw std::invalid_argument(
        "Tensors in assignment operator do not have the same size");
  }
  
  std::copy(rhs.data(), rhs.data() + rhs.n(), s_);

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

Tensor outer(const Vector & a, const Vector & b)
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


} // namespace neml
