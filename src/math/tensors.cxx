#include "tensors.h"

#include "nemlmath.h"

#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace neml {

Tensor::Tensor(std::size_t n) :
    n_(n)
{
  s_ = new double [n_];
}

Tensor::Tensor(const Tensor & other) :
    n_(other.n())
{
  s_ = new double [n_];
  std::copy(other.data(), other.data() + other.n(), s_);
}

Tensor::Tensor(const std::vector<double> flat) : 
  n_(flat.size())
{
  s_ = new double [n_];
  std::copy(flat.begin(), flat.end(), s_);
}

Tensor::Tensor(const double * const flat, size_t n) :
    n_(n)
{
  s_ = new double[n_];
  std::copy(flat, flat + n_, s_);
}

Tensor::~Tensor()
{
  delete [] s_;
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

Tensor & Tensor::operator=(const Tensor && rhs) {
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

Vector::Vector(const double * const v) :
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

RankTwo::RankTwo(const double * const v) :
    Tensor(v, 9)
{
}

RankTwo::RankTwo(const std::vector<const std::vector<double>> A) :
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

} // namespace neml
