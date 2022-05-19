#ifndef TENSORS_H
#define TENSORS_H

#include <vector>
#include <iostream>
#include <map>
#include <cmath>

#include "../windows.h"

namespace neml {

// Forward declarations
class Vector;
class RankTwo;
class Symmetric;
class Skew;
class RankFour;
class SymSymR4;
class SymSkewR4;
class SkewSymR4;
class SymSymSymR6;

class NEML_EXPORT Tensor {
 public:
  Tensor(std::size_t n);
  Tensor(const Tensor & other);
  Tensor(Tensor && other);
  Tensor(const std::vector<double> flat);
  Tensor(double * flat, size_t n);
  Tensor(const double * flat, size_t n);
  virtual ~Tensor();

  /// Do I own my own data?
  bool istore() const {return istore_;};

  Tensor & operator=(const Tensor & rhs);
  Tensor & operator=(Tensor && rhs);

  const double * data() const {return s_;};
  double * s() {return s_;};
  std::size_t n() const {return n_;};

  void copy_data(const double * const indata);

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
NEML_EXPORT bool operator==(const Tensor & a, const Tensor & b);
NEML_EXPORT bool operator!=(const Tensor & a, const Tensor & b);

class NEML_EXPORT Vector: public Tensor {
 public:
  Vector();
  Vector(const std::vector<double> v);
  Vector(double * v);
  Vector(const double * v);

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
NEML_EXPORT Vector operator*(double s, const Vector & v);
NEML_EXPORT Vector operator*(const Vector & v, double s);
NEML_EXPORT Vector operator/(const Vector & v, double s);

// Binary operators with vectors
NEML_EXPORT Vector operator+(const Vector & a, const Vector & b);
NEML_EXPORT Vector operator-(const Vector & a, const Vector & b);

/// io for vectors
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const Vector & v);

// Produces tensors
NEML_EXPORT RankTwo outer(const Vector & a, const Vector & b);

/// Full Rank 2 tensor
class NEML_EXPORT RankTwo: public Tensor {
 public:
  RankTwo();
  RankTwo(const std::vector<double> v);
  RankTwo(double * v);
  RankTwo(const double * v);
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

  static RankTwo id() { return
    RankTwo(std::vector<double>({1.0,0,0,0,1,0,0,0,1}));};

  double norm() const;

  double contract(const RankTwo & other) const;
  double contract(const Symmetric & other) const;
  double contract(const Skew & other) const;
};

// Binary operators with scalars
NEML_EXPORT RankTwo operator*(double s, const RankTwo & v);
NEML_EXPORT RankTwo operator*(const RankTwo & v, double s);
NEML_EXPORT RankTwo operator/(const RankTwo & v, double s);

// Various forms of addition
NEML_EXPORT RankTwo operator+(const RankTwo & a, const RankTwo & b);
NEML_EXPORT RankTwo operator-(const RankTwo & a, const RankTwo & b);
NEML_EXPORT RankTwo operator+(const RankTwo & a, const Symmetric & b);
NEML_EXPORT RankTwo operator-(const RankTwo & a, const Symmetric & b);
NEML_EXPORT RankTwo operator+(const Symmetric & a, const RankTwo & b);
NEML_EXPORT RankTwo operator-(const Symmetric & a, const RankTwo & b);
NEML_EXPORT RankTwo operator+(const RankTwo & a, const Skew & b);
NEML_EXPORT RankTwo operator-(const RankTwo & a, const Skew & b);
NEML_EXPORT RankTwo operator+(const Skew & a, const RankTwo & b);
NEML_EXPORT RankTwo operator-(const Skew & a, const RankTwo & b);

// Various forms of multiplication
NEML_EXPORT Vector operator*(const RankTwo & a, const Vector & b);
NEML_EXPORT Vector operator*(const Vector & a, const RankTwo & b);
NEML_EXPORT RankTwo operator*(const RankTwo & a, const RankTwo & b);
NEML_EXPORT RankTwo operator*(const RankTwo & a, const Symmetric & b);
NEML_EXPORT RankTwo operator*(const Symmetric & a, const RankTwo & b);
NEML_EXPORT RankTwo operator*(const RankTwo & a, const Skew & b);
NEML_EXPORT RankTwo operator*(const Skew & a, const RankTwo & b);
NEML_EXPORT RankTwo operator*(const Skew & a, const Symmetric & b);
NEML_EXPORT RankTwo operator*(const Symmetric & a, const Skew & b);

/// io for RankTwo tensors
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const RankTwo & v);

/// Symmetric Mandel rank 2 tensor
class NEML_EXPORT Symmetric: public Tensor {
 public:
  Symmetric();
  Symmetric(const std::vector<double> v);
  Symmetric(double * v);
  Symmetric(const double * v);
  /// I guess take them seriously and symmetrize it
  Symmetric(const RankTwo & other);

  RankTwo to_full() const;

  Symmetric opposite() const;
  Symmetric operator-() const;

  Symmetric & operator+=(const Symmetric & other);
  Symmetric & operator-=(const Symmetric & other);

  static Symmetric id() { return
    Symmetric(std::vector<double>({1,1,1,0,0,0}));};
  static Symmetric zero() { return
    Symmetric(std::vector<double>({0,0,0,0,0,0}));};
  Symmetric inverse() const;
  Symmetric transpose() const;

  double trace() const;
  Symmetric dev() const;
  double norm() const;

  // Various multiplications
  Vector dot(const Vector & other) const;
  Symmetric dot(const Symmetric & other) const;
  RankTwo dot(const RankTwo & other) const;
  RankTwo dot(const Skew & other) const;

  double contract(const RankTwo & other) const;
  double contract(const Symmetric & other) const;
  double contract(const Skew & other) const;

  double & operator()(size_t i);
  const double & operator()(size_t i) const;
};

// Binary operators with scalars
NEML_EXPORT Symmetric operator*(double s, const Symmetric & v);
NEML_EXPORT Symmetric operator*(const Symmetric & v, double s);
NEML_EXPORT Symmetric operator/(const Symmetric & v, double s);

// Various forms of addition
NEML_EXPORT Symmetric operator+(const Symmetric & a, const Symmetric & b);
NEML_EXPORT Symmetric operator-(const Symmetric & a, const Symmetric & b);

// Various forms of multiplication
NEML_EXPORT Vector operator*(const Symmetric & a, const Vector & b);
NEML_EXPORT Vector operator*(const Vector & a, const Symmetric & b);
NEML_EXPORT Symmetric operator*(const Symmetric & a, const Symmetric & b);

/// io for symmetric tensors
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const Symmetric & v);

class NEML_EXPORT Skew: public Tensor {
 public:
  Skew();
  Skew(const std::vector<double> v);
  Skew(double * v);
  Skew(const double * v);
  /// Skew a general tensor
  Skew(const RankTwo & other);

  static Skew zero() { return
    Skew(std::vector<double>({0,0,0}));};

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

  double contract(const RankTwo & other) const;
  double contract(const Symmetric & other) const;
  double contract(const Skew & other) const;
};

// Binary operators with scalars
NEML_EXPORT Skew operator*(double s, const Skew & v);
NEML_EXPORT Skew operator*(const Skew & v, double s);
NEML_EXPORT Skew operator/(const Skew & v, double s);

// Various forms of addition
NEML_EXPORT Skew operator+(const Skew & a, const Skew & b);
NEML_EXPORT Skew operator-(const Skew & a, const Skew & b);

// Various forms of multiplication
NEML_EXPORT Vector operator*(const Skew & a, const Vector & b);
NEML_EXPORT Vector operator*(const Vector & a, const Skew & b);
NEML_EXPORT Skew operator*(const Skew & a, const Skew & b);

/// io for skew tensors
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const Skew & v);

class NEML_EXPORT RankFour: public Tensor {
 public:
  RankFour();
  RankFour(const std::vector<double> v);
  RankFour(const std::vector<std::vector<std::vector<std::vector<double>>>> A);
  RankFour(double * v);
  RankFour(const double * v);

  RankFour opposite() const;
  RankFour operator-() const;

  RankFour & operator+=(const RankFour & other);
  RankFour & operator-=(const RankFour & other);

  double & operator()(size_t i, size_t j, size_t k, size_t l);
  const double & operator()(size_t i, size_t j, size_t k, size_t l) const;

  SymSymR4 to_sym() const;
  SymSkewR4 to_symskew() const;
  SkewSymR4 to_skewsym() const;

  // Various multiplications
  RankFour dot(const RankFour & other) const;
  RankFour dot(const SymSymR4 & other) const;
  RankFour dot(const SymSkewR4 & other) const;
  RankFour dot(const SkewSymR4 & other) const;

  RankTwo dot(const RankTwo & other) const;
  RankTwo dot(const Symmetric & other) const;
  RankTwo dot(const Skew & other) const;
};

// Binary operators with scalars
NEML_EXPORT RankFour operator*(double s, const RankFour & v);
NEML_EXPORT RankFour operator*(const RankFour & v, double s);
NEML_EXPORT RankFour operator/(const RankFour & v, double s);

// Various forms of addition
NEML_EXPORT RankFour operator+(const RankFour & a, const RankFour & b);
NEML_EXPORT RankFour operator-(const RankFour & a, const RankFour & b);

// Various forms of multiplication
NEML_EXPORT RankFour operator*(const RankFour & a, const RankFour & b);
NEML_EXPORT RankFour operator*(const RankFour & a, const SymSymR4 & b);
NEML_EXPORT RankFour operator*(const RankFour & a, const SymSkewR4 & b);
NEML_EXPORT RankFour operator*(const RankFour & a, const SkewSymR4 & b);

NEML_EXPORT RankTwo operator*(const RankFour & a, const RankTwo & b);
NEML_EXPORT RankTwo operator*(const RankFour & a, const Symmetric & b);
NEML_EXPORT RankTwo operator*(const RankFour & a, const Skew & b);

/// io for SymSymR4 tensors
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const RankFour & v);

class NEML_EXPORT SymSymR4: public Tensor {
 public:
  SymSymR4();
  SymSymR4(const std::vector<double> v);
  SymSymR4(const std::vector<std::vector<double>> A);
  SymSymR4(double * v);
  SymSymR4(const double * v);

  SymSymR4 opposite() const;
  SymSymR4 operator-() const;

  SymSymR4 & operator+=(const SymSymR4 & other);
  SymSymR4 & operator-=(const SymSymR4 & other);

  RankFour to_full() const;

  double & operator()(size_t i, size_t j);
  const double & operator()(size_t i, size_t j) const;

  SymSymR4 inverse() const;
  SymSymR4 transpose() const;

  // Various multiplication
  RankFour dot(const RankFour & other) const;
  SymSymR4 dot(const SymSymR4 & other) const;
  RankFour dot(const SymSkewR4 & other) const;
  RankFour dot(const SkewSymR4 & other) const;

  RankTwo dot(const RankTwo & other) const;
  RankTwo dot(const Skew & other) const;
  Symmetric dot(const Symmetric & other) const;

  static SymSymR4 id() {return SymSymR4(std::vector<std::vector<double>>({
              {1.0,0.0,0.0,0.0,0.0,0.0},
              {0.0,1.0,0.0,0.0,0.0,0.0},
              {0.0,0.0,1.0,0.0,0.0,0.0},
              {0.0,0.0,0.0,1.0,0.0,0.0},
              {0.0,0.0,0.0,0.0,1.0,0.0},
              {0.0,0.0,0.0,0.0,0.0,1.0}
              }
              ));};
  static SymSymR4 id_dev() {return
    SymSymR4(std::vector<std::vector<double>>({
              {1.0-1.0/3.0,-1.0/3.0,-1.0/3.0,0.0,0.0,0.0},
              {-1.0/3.0,1.0-1.0/3.0,-1.0/3.0,0.0,0.0,0.0},
              {-1.0/3.0,-1.0/3.0,1.0-1.0/3.0,0.0,0.0,0.0},
              {0.0,0.0,0.0,1.0,0.0,0.0},
              {0.0,0.0,0.0,0.0,1.0,0.0},
              {0.0,0.0,0.0,0.0,0.0,1.0}
              }
              ));};
  static SymSymR4 zero() {return
    SymSymR4(std::vector<std::vector<double>>({
            {0,0,0,0,0,0},
            {0,0,0,0,0,0},
            {0,0,0,0,0,0},
            {0,0,0,0,0,0},
            {0,0,0,0,0,0},
            {0,0,0,0,0,0}}));};
};

// Binary operators with scalars
NEML_EXPORT SymSymR4 operator*(double s, const SymSymR4 & v);
NEML_EXPORT SymSymR4 operator*(const SymSymR4 & v, double s);
NEML_EXPORT SymSymR4 operator/(const SymSymR4 & v, double s);

// Various forms of addition
NEML_EXPORT SymSymR4 operator+(const SymSymR4 & a, const SymSymR4 & b);
NEML_EXPORT SymSymR4 operator-(const SymSymR4 & a, const SymSymR4 & b);

// Various forms of multiplication
NEML_EXPORT RankFour operator*(const SymSymR4 & a, const RankFour & b);
NEML_EXPORT SymSymR4 operator*(const SymSymR4 & a, const SymSymR4 & b);
NEML_EXPORT RankFour operator*(const SymSymR4 & a, const SymSkewR4 & b);
NEML_EXPORT RankFour operator*(const SymSymR4 & a, const SkewSymR4 & b);

NEML_EXPORT RankTwo operator*(const SymSymR4 & a, const RankTwo & b);
NEML_EXPORT RankTwo operator*(const SymSymR4 & a, const Skew & b);
NEML_EXPORT Symmetric operator*(const SymSymR4 & a, const Symmetric & b);

/// D1_ij D2_kl
NEML_EXPORT SymSymR4 douter(const Symmetric & a, const Symmetric & b);

/// io for SymSymR4 tensors
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const SymSymR4 & v);

class NEML_EXPORT SymSkewR4: public Tensor {
 public:
  SymSkewR4();
  SymSkewR4(const std::vector<double> v);
  SymSkewR4(const std::vector<std::vector<double>> A);
  SymSkewR4(double * v);
  SymSkewR4(const double * v);

  SymSkewR4 opposite() const;
  SymSkewR4 operator-() const;

  SymSkewR4 & operator+=(const SymSkewR4 & other);
  SymSkewR4 & operator-=(const SymSkewR4 & other);

  RankFour to_full() const;

  double & operator()(size_t i, size_t j);
  const double & operator()(size_t i, size_t j) const;

  // Various multiplication
  RankFour dot(const RankFour & other) const;
  RankFour dot(const SymSymR4 & other) const;
  RankFour dot(const SymSkewR4 & other) const;
  RankFour dot(const SkewSymR4 & other) const;

  RankTwo dot(const RankTwo & other) const;
  RankTwo dot(const Skew & other) const;
  RankTwo dot(const Symmetric & other) const;
};

// Binary operators with scalars
NEML_EXPORT SymSkewR4 operator*(double s, const SymSkewR4 & v);
NEML_EXPORT SymSkewR4 operator*(const SymSkewR4 & v, double s);
NEML_EXPORT SymSkewR4 operator/(const SymSkewR4 & v, double s);

// Various forms of addition
NEML_EXPORT SymSkewR4 operator+(const SymSkewR4 & a, const SymSkewR4 & b);
NEML_EXPORT SymSkewR4 operator-(const SymSkewR4 & a, const SymSkewR4 & b);

// Various forms of multiplication
NEML_EXPORT RankFour operator*(const SymSkewR4 & a, const RankFour & b);
NEML_EXPORT RankFour operator*(const SymSkewR4 & a, const SymSymR4 & b);
NEML_EXPORT RankFour operator*(const SymSkewR4 & a, const SymSkewR4 & b);
NEML_EXPORT RankFour operator*(const SymSkewR4 & a, const SkewSymR4 & b);

NEML_EXPORT RankTwo operator*(const SymSkewR4 & a, const RankTwo & b);
NEML_EXPORT RankTwo operator*(const SymSkewR4 & a, const Skew & b);
NEML_EXPORT RankTwo operator*(const SymSkewR4 & a, const Symmetric & b);

/// io for SymSkewR4 tensors
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const SymSkewR4 & v);

class NEML_EXPORT SkewSymR4: public Tensor {
 public:
  SkewSymR4();
  SkewSymR4(const std::vector<double> v);
  SkewSymR4(const std::vector<std::vector<double>> A);
  SkewSymR4(double * v);
  SkewSymR4(const double * v);

  SkewSymR4 opposite() const;
  SkewSymR4 operator-() const;

  SkewSymR4 & operator+=(const SkewSymR4 & other);
  SkewSymR4 & operator-=(const SkewSymR4 & other);

  RankFour to_full() const;

  double & operator()(size_t i, size_t j);
  const double & operator()(size_t i, size_t j) const;

  // Various multiplication
  RankFour dot(const RankFour & other) const;
  RankFour dot(const SymSymR4 & other) const;
  RankFour dot(const SkewSymR4 & other) const;
  RankFour dot(const SymSkewR4 & other) const;

  RankTwo dot(const RankTwo & other) const;
  RankTwo dot(const Skew & other) const;
  RankTwo dot(const Symmetric & other) const;
};

// Binary operators with scalars
NEML_EXPORT SkewSymR4 operator*(double s, const SkewSymR4 & v);
NEML_EXPORT SkewSymR4 operator*(const SkewSymR4 & v, double s);
NEML_EXPORT SkewSymR4 operator/(const SkewSymR4 & v, double s);

// Various forms of addition
NEML_EXPORT SkewSymR4 operator+(const SkewSymR4 & a, const SkewSymR4 & b);
NEML_EXPORT SkewSymR4 operator-(const SkewSymR4 & a, const SkewSymR4 & b);

// Various forms of multiplication
NEML_EXPORT RankFour operator*(const SkewSymR4 & a, const RankFour & b);
NEML_EXPORT RankFour operator*(const SkewSymR4 & a, const SymSymR4 & b);
NEML_EXPORT RankFour operator*(const SkewSymR4 & a, const SkewSymR4 & b);
NEML_EXPORT RankFour operator*(const SkewSymR4 & a, const SymSkewR4 & b);

NEML_EXPORT RankTwo operator*(const SkewSymR4 & a, const RankTwo & b);
NEML_EXPORT RankTwo operator*(const SkewSymR4 & a, const Skew & b);
NEML_EXPORT RankTwo operator*(const SkewSymR4 & a, const Symmetric & b);

// ij kl
NEML_EXPORT SkewSymR4 douter(const Skew & a, const Symmetric & b);

/// io for SkewSymR4 tensors
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const SkewSymR4 & v);


// Helper CP methods
/// Skmab*Wml - Wkm*Smlab
NEML_EXPORT SymSymR4 SymSymR4Skew_SkewSymR4SymR4(const SymSymR4 & S, const Skew & W);

/// Dkm*Smlab - Skmab*Dml
NEML_EXPORT SymSymR4 SymSkewR4Sym_SkewSymR4SymR4(const SkewSymR4 & S, const Symmetric & D);

/// Specialty operator for Skew part C_ijkb e_ka - C_ijal e_bl
NEML_EXPORT SymSkewR4 SpecialSymSymR4Sym(const SymSymR4 & S, const Symmetric & D);

/// Hopefully the only outer product I'll need...
NEML_EXPORT SymSymSymR6 outer_product_k(const SymSymR4 & A, const Symmetric & B);

/// Rank six where all adjacent indices have minor symmetry
class NEML_EXPORT SymSymSymR6: public Tensor {
 public:
  SymSymSymR6();
  SymSymSymR6(const std::vector<double> v);
  SymSymSymR6(const std::vector<std::vector<std::vector<double>>> A);
  SymSymSymR6(double * v);
  SymSymSymR6(const double * v);

  SymSymSymR6 opposite() const;
  SymSymSymR6 operator-() const;

  SymSymSymR6 & operator+=(const SymSymSymR6 & other);
  SymSymSymR6 & operator-=(const SymSymSymR6 & other);

  double & operator()(size_t i, size_t j, size_t k);
  const double & operator()(size_t i, size_t j, size_t k) const;

  SymSymR4 dot_i(const Symmetric & other) const;
  SymSymR4 dot_j(const Symmetric & other) const;
  SymSymR4 dot_k(const Symmetric & other) const;

  SymSymSymR6 middle_dot_after(const SymSymR4 & other) const;
  SymSymSymR6 middle_dot_before(const SymSymR4 & other) const;
};

// Binary operators with scalars
NEML_EXPORT SymSymSymR6 operator*(double s, const SymSymSymR6 & v);
NEML_EXPORT SymSymSymR6 operator*(const SymSymSymR6 & v, double s);
NEML_EXPORT SymSymSymR6 operator/(const SymSymSymR6 & v, double s);

// Various forms of addition
NEML_EXPORT SymSymSymR6 operator+(const SymSymSymR6 & a, const SymSymSymR6 & b);
NEML_EXPORT SymSymSymR6 operator-(const SymSymSymR6 & a, const SymSymSymR6 & b);

} // namespace neml

#endif
