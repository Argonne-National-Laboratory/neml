#ifndef CRYSTALLOGRAPHY_H
#define CRYSTALLOGRAPHY_H

#include "../objects.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include "../windows.h"

#include <vector>
#include <string>
#include <memory>

namespace neml {

std::vector<Orientation> symmetry_rotations(std::string sclass);

class SymmetryGroup: public NEMLObject {
 public:
  /// Initialize with the Hermann-Mauguin notation as a string
  SymmetryGroup(std::string sclass);
  /// Destructor
  virtual ~SymmetryGroup();

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Quaternion symmetry operators
  const std::vector<Orientation> & ops() const;

  /// Number of symmetry operators
  size_t nops() const;

  /// Find the minimum misorientation transformation between a and b
  Orientation misorientation(const Orientation & a, const Orientation & b) const;

  /// Find the distance between a and b with the usual metric
  double distance(const Orientation & a, const Orientation & b) const;

  /// Find the orientation which, under symmetry operations, is closest to b
  void closest(const Orientation & a, const Orientation & b,
               Orientation & close, double & dist) const;

 private:
  void opdist_(const double * const q1, const double * const q2,
               double * const res) const;

 private:
  const std::vector<Orientation> ops_;
  std::vector<double> raw_;
};

static Register<SymmetryGroup> regSymmetryGroup;

class Lattice: public NEMLObject {
 public:
  /// Initialize with the three lattice vectors, the symmetry group and
  /// (optionally) a initial list of slip systems
  Lattice(Vector a1, Vector a2, Vector a3, std::shared_ptr<SymmetryGroup> symmetry,
          list_systems isystems = {});
  /// Destructor
  virtual ~Lattice();

  /// First lattice vector
  const Vector & a1() {return a1_;};
  /// Second lattice vector
  const Vector & a2() {return a2_;};
  /// Third lattice vector
  const Vector & a3() {return a3_;};
  /// First reciprocal vector
  const Vector & b1() {return b1_;};
  /// Second reciprocal vector
  const Vector & b2() {return b2_;};
  /// Third reciprocal vector
  const Vector & b3() {return b3_;};

  /// Return the list of burgers vectors
  const std::vector<std::vector<Vector>> & burgers_vectors() {return burgers_vectors_;};
  /// Return the list of normalized slip directions
  const std::vector<std::vector<Vector>> & slip_directions() {return slip_directions_;};
  /// Return the last of normalize slip normals
  const std::vector<std::vector<Vector>> & slip_planes() {return slip_planes_;};

  /// Convert Miller directions to cartesian vectors
  Vector miller2cart_direction(std::vector<int> m);
  /// Convert Miller planes to cartesian normal vectors
  Vector miller2cart_plane(std::vector<int> m);

  /// Find all sets of equivalent vectors (+/- different)
  std::vector<Vector> equivalent_vectors(Vector v);
  /// Find all all sets of equivalent vectors (+/- the same)
  std::vector<Vector> equivalent_vectors_bidirectional(Vector v);

  /// Add a slip system given the Miller direction and plane
  void add_slip_system(std::vector<int> d, std::vector<int> p);

  /// Number of groups of slip systems
  size_t ngroup() const;
  /// Number of slip systems in group g
  size_t nslip(size_t g) const;
  /// Flat index of slip group g, system i
  size_t flat(size_t g, size_t i) const;

  /// Return the sym(d x n) tensor for group g, system i, rotated with Q
  Symmetric M(size_t g, size_t i, const Orientation & Q);
  /// Return the skew(d x n) tensor for group g, system i, rotated with Q
  Skew N(size_t g, size_t i, const Orientation & Q);

  /// Calculate the resolved shear stress on group g, system i, rotated with Q
  /// given the stress
  double shear(size_t g, size_t i, const Orientation & Q, const Symmetric &
               stress);
  /// Calculate the derivative of the resolved shear stress on group g,
  /// system i, rotated with Q, given the stress
  Symmetric d_shear(size_t g, size_t i, const Orientation & Q, const Symmetric &
                    stress);

  /// Access the symmetry operations
  const std::shared_ptr<SymmetryGroup> symmetry();

 private:
  void make_reciprocal_lattice_();
  static void assert_miller_(std::vector<int> m);

  void cache_rot_(const Orientation & Q);

 private:
  Vector a1_, a2_, a3_, b1_, b2_, b3_;
  const std::shared_ptr<SymmetryGroup> symmetry_;

  std::vector<std::vector<Vector>> burgers_vectors_;
  std::vector<std::vector<Vector>> slip_directions_;
  std::vector<std::vector<Vector>> slip_planes_;

  std::vector<size_t> offsets_;

  // Used for caching common asks
  bool setup_;
  size_t hash_;
  std::vector<std::vector<Symmetric>> Ms_;
  std::vector<std::vector<Skew>> Ns_;
};

class CubicLattice: public Lattice {
 public:
  /// Specialized Lattice for cubic systems, initialize with the lattice
  /// parameter
  CubicLattice(double a, list_systems isystems = {});

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();
};

static Register<CubicLattice> regCubicLattice;

} // namespace neml

#endif // CRYSTALLOGRAPHY_H
