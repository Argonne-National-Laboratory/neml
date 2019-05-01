#ifndef TENSORS_H
#define TENSORS_H

#include <vector>
#include <iostream>

namespace neml {

class Vector;
class RankTwo;

class Tensor {
 public:
  Tensor(std::size_t n);
  Tensor(const Tensor & other);
  Tensor(const std::vector<double> flat);
  Tensor(const double * const flat, size_t n);
  virtual ~Tensor();

  Tensor & operator=(const Tensor & rhs);
  Tensor & operator=(const Tensor && rhs);

  const double * data() const {return s_;};
  double * s() {return s_;};
  std::size_t n() const {return n_;};

  Tensor & operator*=(double s);
  Tensor & operator/=(double s);

 protected:
  /// Helper to add data
  void add_(const Tensor & other);
  
  /// Helper to negate
  void negate_();

 protected:
  double * s_;
  const std::size_t n_;
};

/// Dangerous but useful
bool operator==(const Tensor & a, const Tensor & b);
bool operator!=(const Tensor & a, const Tensor & b);

class Vector: public Tensor {
 public:
  Vector();
  Vector(const std::vector<double> v);
  Vector(const double * const v);

  Vector opposite() const;
  Vector operator-() const;

  Vector & operator+=(const Vector & other);
  Vector & operator-=(const Vector & other);

  double & operator()(size_t i);
  const double & operator()(size_t i) const;

  double dot(const Vector & rhs) const;
  RankTwo outer(const Vector & o) const;

  double norm() const;
  Vector cross(const Vector & other) const;
  Vector & normalize();
};

// Binary operators with scalars
Vector operator*(double s, const Vector & v);
Vector operator*(const Vector & v, double s);
Vector operator/(const Vector & v, double s);

// Binary operators with vectors
Vector operator+(const Vector & a, const Vector & b);
Vector operator-(const Vector & a, const Vector & b);

/// io for vectors
std::ostream & operator<<(std::ostream & os, const Vector & v);

// Produces tensors
Tensor outer(const Vector & a, const Vector & b);

/// Full Rank 2 tensor
class RankTwo: public Tensor {
 public:
  RankTwo();
  RankTwo(const std::vector<double> v);
  RankTwo(const double * const v);
  RankTwo(const std::vector<const std::vector<double>> A);

  RankTwo opposite() const;
  RankTwo operator-() const;

  RankTwo & operator+=(const RankTwo & other);
  RankTwo & operator-=(const RankTwo & other);

  double & operator()(size_t i, size_t j);
  const double & operator()(size_t i, size_t j) const;
};

// Binary operators with scalars
RankTwo operator*(double s, const RankTwo & v);
RankTwo operator*(const RankTwo & v, double s);
RankTwo operator/(const RankTwo & v, double s);

// Binary operators with vectors
RankTwo operator+(const RankTwo & a, const RankTwo & b);
RankTwo operator-(const RankTwo & a, const RankTwo & b);

/// io for tensors
std::ostream & operator<<(std::ostream & os, const RankTwo & v);

} // namespace neml

#endif
