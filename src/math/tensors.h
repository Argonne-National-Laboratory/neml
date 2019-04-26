#ifndef TENSORS_H
#define TENSORS_H

#include <vector>
#include <iostream>

namespace neml {

class Tensor {
 public:
  Tensor(std::size_t n);
  Tensor(const Tensor & other);
  virtual ~Tensor();

  Tensor & operator=(const Tensor & rhs);
  Tensor & operator=(const Tensor && rhs);

  const double * data() const {return s_;};
  double * s() {return s_;};
  std::size_t n() const {return n_;};

 protected:
  double * s_;
  const std::size_t n_;
};

class Vector: public Tensor {
 public:
  Vector();
  Vector(const std::vector<double> v);
  Vector(const double * const v);

  double dot(const Vector & rhs) const;
  double norm() const;

  Vector opposite() const;
  Vector operator-() const;

  Vector & operator*=(double s);
  Vector & operator/=(double s);

  Vector & operator+=(const Vector & other);
  Vector & operator-=(const Vector & other);

  Vector cross(const Vector & other) const;

  Vector & normalize();
};

// Binary operators
Vector operator*(double s, const Vector & v);
Vector operator*(const Vector & v, double s);
Vector operator/(double s, const Vector & v);
Vector operator/(const Vector & v, double s);

Vector operator+(const Vector & a, const Vector & b);
Vector operator-(const Vector & a, const Vector & b);

/// Dangerous but useful
bool operator==(const Vector & a, const Vector & b);
bool operator!=(const Vector & a, const Vector & b);

/// io for vectors
std::ostream & operator<<(std::ostream & os, const Vector & v);

} // namespace neml

#endif
