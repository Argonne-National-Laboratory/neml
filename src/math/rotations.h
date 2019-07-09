#ifndef ROTATIONS_H
#define ROTATIONS_H

#include "tensors.h"

#include "../objects.h"

#include <complex>
#include <vector>
#include <string>

#define CINDEX(i,j,n) (j + i * n)

namespace neml {

/// A generic quaternion, stored as [s v1 v2 v3]
class Quaternion {
 public:
  Quaternion();
  Quaternion(const std::vector<double> v);
  Quaternion(double * v);
  Quaternion(const Quaternion & other);
  Quaternion(const Quaternion && other);

  virtual ~Quaternion();
  
  /// Copy
  Quaternion & operator=(const Quaternion & rhs);
  /// Move
  Quaternion & operator=(const Quaternion && rhs);
  /// Do you store your own data?
  bool store() const;

  /// Raw quaternion
  const double * quat() const;
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
Quaternion operator*(double s, const Quaternion & q);
Quaternion operator*(const Quaternion & q, double s);
Quaternion operator*(const Quaternion & lhs, const Quaternion & rhs);
Quaternion operator/(const Quaternion & q, double s);
Quaternion operator/(const Quaternion & lhs, const Quaternion & rhs);

class Orientation: public Quaternion, public NEMLObject {
 public:
  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  // Creation functions
  /// Create from a Rodrigues vector
  static Orientation createRodrigues(const double * const r);
  void setRodrigues(const double * const r);

  /// Create from a rotation matrix
  static Orientation createMatrix(const double * const M);
  void setMatrix(const double * const M);

  /// Create from an axis-angle representation
  static Orientation createAxisAngle(const double * const n, double a,
                                     std::string angles = "radians");
  void setAxisAngle(const double * const n, double a, 
                    std::string angles = "radians");

  /// Create from various Euler angles
  static Orientation createEulerAngles(double a, double b, double c, 
                                       std::string angles = "radians",
                                       std::string convention = "kocks");
  void setEulerAngles(double a, double b, double c, 
                      std::string angles = "radians",
                      std::string convention = "kocks");

  /// Create from the Hopf coordinates
  static Orientation createHopf(double psi, double theta, double phi,
                                std::string angles = "radians");
  void setHopf(double psi, double theta, double phi, 
               std::string angles = "radians");

  /// Create from hyperspherical coordinates
  static Orientation createHyperspherical(double a1, double a2, double a3,
                                          std::string angles = "radians");
  void setHyperspherical(double a1, double a2, double a3, 
                         std::string angles = "radians");
  
  // Actual constructors
  Orientation();
  Orientation(double * v);
  Orientation(const std::vector<double> v);
  Orientation(const Quaternion & other);

  // Various conversions
  void to_euler(double & a, double & b, double & c, 
                std::string angles = "radians",
                std::string convention = "kocks") const;
  void to_axis_angle(double * const n, double & a,
                     std::string angles = "radians") const;
  void to_matrix(double * const M) const;
  RankTwo to_tensor() const;
  void to_rodrigues(double * const v) const;
  void to_hopf(double & alpha, double & beta, double & gamma,
               std::string angles = "radians") const;
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

  /// Rotate various tensors
  Vector apply(const Vector & a) const;
  RankTwo apply(const RankTwo & a) const;
  Symmetric apply(const Symmetric & a) const;
  Skew apply(const Skew & a) const;
  RankFour apply(const RankFour & a) const;
  SymSym apply(const SymSym & a) const;

 private:
  void normalize_();
  static void to_kocks_(double a, double b, double c, double & oa, double & ob,
                        double & oc, std::string type);
  static void from_kocks_(double a, double b, double c, double & oa, double & ob,
                          double & oc, std::string type);
  static void kocks_to_matrix_(double a, double b, double c, double * const M);
};

static Register<Orientation> regOrientation;

// Binary operators
Orientation operator*(const Orientation & lhs, const Orientation & rhs);
Orientation operator/(const Orientation & lhs, const Orientation & rhs);

/// Generate n random orientations
//    This algorithm comes from LaValle, 2006 who I believe grabbed it
//    from Shoemake, 1992
std::vector<Orientation> random_orientations(int n);

/// Exponential map of a skew tensor in my convention
Orientation wexp(const Skew & w);

/// Inverse exponential map of a quaternion to a skew tensor in my convention
Skew wlog(const Orientation & q);

} // namespace neml

#endif // ROTATIONS_H
