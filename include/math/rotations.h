#ifndef ROTATIONS_H
#define ROTATIONS_H

#include "tensors.h"
#include "nemlmath.h"

#include "../objects.h"

#include <complex>
#include <vector>
#include <string>

#include "../windows.h"

namespace neml {

/// A generic quaternion, stored as [s v1 v2 v3]
class NEML_EXPORT Quaternion {
 public:
  /// Default constructor (manage own memory)
  Quaternion();
  /// Construct from vector (manage own memory)
  Quaternion(const std::vector<double> v);
  /// Construct from raw pointer (don't manage memory)
  Quaternion(double * v);
  /// Copy constructor
  Quaternion(const Quaternion & other);
  /// Move constructor
  Quaternion(const Quaternion && other);

  virtual ~Quaternion();

  /// Copy
  Quaternion & operator=(const Quaternion & rhs);
  /// Move
  Quaternion & operator=(const Quaternion && rhs);
  /// Do you store your own data?
  bool store() const;

  /// Raw quaternion as a const pointer
  const double * quat() const;
  /// Raw quaternion as a nonconst pointer
  double * data();

  /// Quaternion norm
  virtual double norm() const;
  /// Opposite
  Quaternion opposite() const;
  /// C++ operator opposite
  Quaternion operator-() const;
  /// Conjugation
  Quaternion conj() const;
  /// Opposite scalar
  Quaternion flip() const;
  /// Inversion
  Quaternion inverse() const;
  /// Exponential map
  Quaternion exp() const;
  /// Inverse exponential map
  Quaternion log() const;
  /// Quaternion multiplication
  Quaternion & operator*=(const Quaternion & rhs);
  /// Scalar multiplication
  Quaternion & operator*=(double scalar);
  /// Quaternion division
  Quaternion & operator/=(const Quaternion & rhs);
  /// Scalar division
  Quaternion & operator/=(double scalar);

  /// Power
  Quaternion pow(double w) const;

  /// Dot product, useful for various distances
  double dot(const Quaternion & other) const;

  /// Hash function for quick comparisons
  size_t hash() const;

  /// Product matrix for quaternion composition
  void to_product_matrix(double * M) const;

 protected:
  void alloc_();

 protected:
  void opposite_(double * const out) const;
  void conj_(double * const out) const;
  void flip_(double * const out) const;
  void inverse_(double * const out) const;
  void exp_(double * const out) const;
  void log_(double * const out) const;
  void multiply_(const double * const b);
  void smultiply_(double s);

  double * quat_;
  bool store_;
};

// Binary operators
/// Scalar multiplication
NEML_EXPORT Quaternion operator*(double s, const Quaternion & q);
/// Scalar multiplication
NEML_EXPORT Quaternion operator*(const Quaternion & q, double s);
/// Composition
NEML_EXPORT Quaternion operator*(const Quaternion & lhs, const Quaternion & rhs);
/// Scalar division
NEML_EXPORT Quaternion operator/(const Quaternion & q, double s);
/// Scalar division
NEML_EXPORT Quaternion operator/(const Quaternion & lhs, const Quaternion & rhs);

/// C++ stream output for quaternions
NEML_EXPORT std::ostream & operator<<(std::ostream & os, const Quaternion & q);

class NEML_EXPORT Orientation: public Quaternion {
 public:
  // Creation functions
  /// Create from a Rodrigues vector
  static Orientation createRodrigues(const double * const r);
  /// Set from an input Rodrigues vector
  void setRodrigues(const double * const r);

  /// Create from a rotation matrix
  static Orientation createMatrix(const double * const M);
  /// Set from an input matrix
  void setMatrix(const double * const M);

  /// Create from an axis-angle representation
  static Orientation createAxisAngle(const double * const n, double a,
                                     std::string angles = "radians");
  /// Set from an input axis/angle pair
  void setAxisAngle(const double * const n, double a,
                    std::string angles = "radians");

  /// Create from various Euler angles
  static Orientation createEulerAngles(double a, double b, double c,
                                       std::string angles = "radians",
                                       std::string convention = "kocks");
  /// Set from inputEuler angles
  void setEulerAngles(double a, double b, double c,
                      std::string angles = "radians",
                      std::string convention = "kocks");

  /// Create from the Hopf coordinates
  static Orientation createHopf(double psi, double theta, double phi,
                                std::string angles = "radians");
  /// Set from input Hopf coordinates
  void setHopf(double psi, double theta, double phi,
               std::string angles = "radians");

  /// Create from hyperspherical coordinates
  static Orientation createHyperspherical(double a1, double a2, double a3,
                                          std::string angles = "radians");
  /// Set from input hyperspherical coordinates
  void setHyperspherical(double a1, double a2, double a3,
                         std::string angles = "radians");

  /// Create from two vectors
  static Orientation createVectors(const Vector & x, const Vector & y);
  /// Set from input two vectors
  void setVectors(const Vector & x, const Vector & y);

  // Actual constructors
  /// Default constructor (defaults to identity, manage own memory)
  Orientation();
  /// Raw pointer constructor (don't manage memory)
  Orientation(double * v);
  /// vector<double> constructor (manage own memory)
  Orientation(const std::vector<double> v);
  /// Copy constructor
  Orientation(const Quaternion & other);

  /// Explicit deepcopy
  Orientation deepcopy() const;

  // Various conversions
  /// Convert to Euler angles
  void to_euler(double & a, double & b, double & c,
                std::string angles = "radians",
                std::string convention = "kocks") const;
  /// Convert to an axis/angle pair
  void to_axis_angle(double * const n, double & a,
                     std::string angles = "radians") const;
  /// Convert to a rotation matrix
  void to_matrix(double * const M) const;
  /// Convert to a rank 2 tensor
  RankTwo to_tensor() const;
  /// Convert to a Rodrigues vector
  void to_rodrigues(double * const v) const;
  /// Convert to Hopf coordinates
  void to_hopf(double & alpha, double & beta, double & gamma,
               std::string angles = "radians") const;
  /// Convert to hyperspherical coordinates
  void to_hyperspherical(double & a1, double & a2, double & a3,
                         std::string angles = "radians") const;

  // Annoyingly the operators must be for the most part redefined
  /// Opposite
  Orientation opposite() const;
  /// C++ operator opposite
  Orientation operator-() const;
  /// Conjugation
  Orientation conj() const;
  /// Opposite scalar
  Orientation flip() const;
  /// Inversion
  Orientation inverse() const;
  /// Orientation multiplication
  Orientation & operator*=(const Orientation & rhs);
  /// Orientation division
  Orientation & operator/=(const Orientation & rhs);

  /// Power
  Orientation pow(double w) const;

  /// Rotate a Vector
  Vector apply(const Vector & a) const;
  /// Rotate a Rank2
  RankTwo apply(const RankTwo & a) const;
  /// Rotate a Symmetric tensor
  Symmetric apply(const Symmetric & a) const;
  /// Rotate a Skew tensor
  Skew apply(const Skew & a) const;
  /// Rotate a RankFour tensor
  RankFour apply(const RankFour & a) const;
  /// Rotate a SymSymR4 rank four tensor
  SymSymR4 apply(const SymSymR4 & a) const;

  /// Geodesic distance
  double distance(const Orientation & other) const;

 private:
  void normalize_();
  static void to_kocks_(double a, double b, double c, double & oa, double & ob,
                        double & oc, std::string type);
  static void from_kocks_(double a, double b, double c, double & oa, double & ob,
                          double & oc, std::string type);
  static void kocks_to_matrix_(double a, double b, double c, double * const M);
};

// Binary operators
/// Compose two rotations
NEML_EXPORT Orientation operator*(const Orientation & lhs, const Orientation & rhs);
/// Compose a rotation with the inverse of a rotation
NEML_EXPORT Orientation operator/(const Orientation & lhs, const Orientation & rhs);

class NEML_EXPORT CrystalOrientation : public NEMLObject, public Orientation { 
 public:
  CrystalOrientation(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
};

static Register<CrystalOrientation> regCrystalOrientation;

std::shared_ptr<CrystalOrientation> zero_orientation();

/// Generate n random orientations
//    This algorithm comes from LaValle, 2006 who I believe grabbed it
//    from Shoemake, 1992
NEML_EXPORT std::vector<CrystalOrientation> random_orientations(int n);

/// Exponential map of a skew tensor in my convention
NEML_EXPORT Orientation wexp(const Skew & w);

/// Inverse exponential map of a quaternion to a skew tensor in my convention
NEML_EXPORT Skew wlog(const Orientation & q);

/// Geodesic distance
NEML_EXPORT double distance(const Orientation & q1, const Orientation & q2);

/// Arbitrary rotation from a to b
NEML_EXPORT Orientation rotate_to(const Vector & a, const Vector & b);

/// Family of rotations from a to b parameterized by an angle
NEML_EXPORT Orientation rotate_to_family(const Vector & a, const Vector & b, double ang);

/// Convert an orientation to a crystal orientation
NEML_EXPORT CrystalOrientation make_crystal_orientation(const Orientation & o);

} // namespace neml

#endif // ROTATIONS_H