#ifndef TENSORS_H
#define TENSORS_H

#include <vector>
#include <iostream>
#include <map>
#include <cmath>

namespace neml {

class Vector;
class RankTwo;
class Symmetric;
class Skew;

class Tensor {
 public:
  Tensor(std::size_t n);
  Tensor(const Tensor & other);
  Tensor(const Tensor && other);
  Tensor(const std::vector<double> flat);
  Tensor(double * flat, size_t n);
  virtual ~Tensor();
  
  /// Do I own my own data?
  bool istore() const {return istore_;};

  Tensor & operator=(const Tensor & rhs);
  Tensor & operator=(Tensor && rhs);

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
  bool istore_;
};

/// Dangerous but useful
bool operator==(const Tensor & a, const Tensor & b);
bool operator!=(const Tensor & a, const Tensor & b);

class Vector: public Tensor {
 public:
  Vector();
  Vector(const std::vector<double> v);
  Vector(double * v);

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
  RankTwo(double * v);
  RankTwo(const std::vector<std::vector<double>> A);
  /// Helper
  RankTwo(const Symmetric & other);
  RankTwo(const Skew & other);

  RankTwo opposite() const;
  RankTwo operator-() const;

  RankTwo & operator+=(const RankTwo & other);
  RankTwo & operator-=(const RankTwo & other);

  RankTwo & operator+=(const Symmetric & other);
  RankTwo & operator-=(const Symmetric & other);

  RankTwo & operator+=(const Skew & other);
  RankTwo & operator-=(const Skew & other);

  double & operator()(size_t i, size_t j);
  const double & operator()(size_t i, size_t j) const;

  RankTwo dot(const RankTwo & other) const;
  RankTwo dot(const Symmetric & other) const;
  RankTwo dot(const Skew & other) const;
  Vector dot(const Vector & other) const;

  RankTwo inverse() const;
  RankTwo transpose() const;
};

// Binary operators with scalars
RankTwo operator*(double s, const RankTwo & v);
RankTwo operator*(const RankTwo & v, double s);
RankTwo operator/(const RankTwo & v, double s);

// Various forms of addition
RankTwo operator+(const RankTwo & a, const RankTwo & b);
RankTwo operator-(const RankTwo & a, const RankTwo & b);
RankTwo operator+(const RankTwo & a, const Symmetric & b);
RankTwo operator-(const RankTwo & a, const Symmetric & b);
RankTwo operator+(const Symmetric & a, const RankTwo & b);
RankTwo operator-(const Symmetric & a, const RankTwo & b);
RankTwo operator+(const RankTwo & a, const Skew & b);
RankTwo operator-(const RankTwo & a, const Skew & b);
RankTwo operator+(const Skew & a, const RankTwo & b);
RankTwo operator-(const Skew & a, const RankTwo & b);

// Various forms of multiplication
Vector operator*(const RankTwo & a, const Vector & b);
Vector operator*(const Vector & a, const RankTwo & b);
RankTwo operator*(const RankTwo & a, const RankTwo & b);
RankTwo operator*(const RankTwo & a, const Symmetric & b);
RankTwo operator*(const Symmetric & a, const RankTwo & b);
RankTwo operator*(const RankTwo & a, const Skew & b);
RankTwo operator*(const Skew & a, const RankTwo & b);
RankTwo operator*(const Skew & a, const Symmetric & b);
RankTwo operator*(const Symmetric & a, const Skew & b);

/// io for RankTwo tensors
std::ostream & operator<<(std::ostream & os, const RankTwo & v);

/// Symmetric Mandel rank 2 tensor
class Symmetric: public Tensor {
 public:
  Symmetric();
  Symmetric(const std::vector<double> v);
  Symmetric(double * v);
  /// I guess take them seriously and symmetrize it
  Symmetric(const RankTwo & other);

  RankTwo to_full() const;

  Symmetric opposite() const;
  Symmetric operator-() const;

  Symmetric & operator+=(const Symmetric & other);
  Symmetric & operator-=(const Symmetric & other);

  Symmetric inverse() const;
  Symmetric transpose() const;

  // Various multiplications
  Vector dot(const Vector & other) const;
  Symmetric dot(const Symmetric & other) const;
  RankTwo dot(const RankTwo & other) const;
  RankTwo dot(const Skew & other) const;
};

// Binary operators with scalars
Symmetric operator*(double s, const Symmetric & v);
Symmetric operator*(const Symmetric & v, double s);
Symmetric operator/(const Symmetric & v, double s);

// Various forms of addition
Symmetric operator+(const Symmetric & a, const Symmetric & b);
Symmetric operator-(const Symmetric & a, const Symmetric & b);

// Various forms of multiplication
Vector operator*(const Symmetric & a, const Vector & b);
Vector operator*(const Vector & a, const Symmetric & b);
Symmetric operator*(const Symmetric & a, const Symmetric & b);

/// io for symmetric tensors
std::ostream & operator<<(std::ostream & os, const Symmetric & v);

class Skew: public Tensor {
 public:
  Skew();
  Skew(const std::vector<double> v);
  Skew(double * v);
  /// Skew a general tensor
  Skew(const RankTwo & other);

  RankTwo to_full() const;

  Skew opposite() const;
  Skew operator-() const;

  Skew & operator+=(const Skew & other);
  Skew & operator-=(const Skew & other);

  Skew transpose() const;

  // Various multiplication
  Vector dot(const Vector & other) const;
  Skew dot(const Skew & other) const;
  RankTwo dot(const RankTwo & other) const;
  RankTwo dot(const Symmetric & other) const;
};

// Binary operators with scalars
Skew operator*(double s, const Skew & v);
Skew operator*(const Skew & v, double s);
Skew operator/(const Skew & v, double s);

// Various forms of addition
Skew operator+(const Skew & a, const Skew & b);
Skew operator-(const Skew & a, const Skew & b);

// Various forms of multiplication
Vector operator*(const Skew & a, const Vector & b);
Vector operator*(const Vector & a, const Skew & b);
Skew operator*(const Skew & a, const Skew & b);

/// io for symmetric tensors
std::ostream & operator<<(std::ostream & os, const Skew & v);

} // namespace neml

#endif
