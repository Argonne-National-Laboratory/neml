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
    Tensor(3)
{
  if (v.size() != 3) {
    throw std::invalid_argument("Input to vector must have size 3!");
  }
  std::copy(v.begin(), v.end(), s_);
}

Vector::Vector(const double * const v) :
    Tensor(3)
{
  std::copy(v, v+3, s_);
}

double Vector::dot(const Vector & rhs) const
{
  double s = 0.0;
  for (std::size_t i = 0; i < 3; i++) {
    s += s_[i] * rhs.data()[i];
  }
  return s;
}

double Vector::norm() const
{
  return sqrt(dot(*this));
}

Vector Vector::opposite() const
{
  double s[3];
  for (int i=0; i<3; i++) {
    s[i] = -s_[i];
  }
  return Vector(s);
}

Vector Vector::operator-() const
{
  return opposite();
}

Vector & Vector::operator+=(const Vector & other)
{
  for (int i=0; i<3; i++) {
    s_[i] += other.data()[i];
  }

  return *this;
}

Vector & Vector::operator-=(const Vector & other)
{
  return this->operator+=(-other);
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

} // namespace neml
