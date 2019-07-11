#ifndef CRYSTALLOGRAPHY_H
#define CRYSTALLOGRAPHY_H

#include "../objects.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include <vector>
#include <string>
#include <memory>

namespace neml {

extern const std::vector<Orientation> tetragonal;
extern const std::vector<Orientation> hexagonal;
extern const std::vector<Orientation> cubic;

std::vector<Orientation> symmetry_rotations(std::string sclass);

class SymmetryGroup: public NEMLObject {
 public:
  SymmetryGroup(std::string sclass);
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
  Lattice(Vector a1, Vector a2, Vector a3, std::shared_ptr<SymmetryGroup> symmetry,
          list_systems isystems = {});
  virtual ~Lattice();

  const Vector & a1() {return a1_;};
  const Vector & a2() {return a2_;};
  const Vector & a3() {return a3_;};
  const Vector & b1() {return b1_;};
  const Vector & b2() {return b2_;};
  const Vector & b3() {return b3_;};

  const std::vector<std::vector<Vector>> & burgers_vectors() {return burgers_vectors_;};
  const std::vector<std::vector<Vector>> & slip_directions() {return slip_directions_;};
  const std::vector<std::vector<Vector>> & slip_planes() {return slip_planes_;};
  
  Vector miller2cart_direction(std::vector<int> m);
  Vector miller2cart_plane(std::vector<int> m);

  std::vector<Vector> equivalent_vectors(Vector v);
  std::vector<Vector> equivalent_vectors_bidirectional(Vector v);

  void add_slip_system(std::vector<int> d, std::vector<int> p);

  size_t ngroup() const;
  size_t nslip(size_t g) const;
  size_t flat(size_t g, size_t i) const;

  Symmetric M(size_t g, size_t i, const Orientation & Q) const;
  Skew N(size_t g, size_t i, const Orientation & Q) const;

  double shear(size_t g, size_t i, const Orientation & Q, const Symmetric &
               stress) const;
  Symmetric d_shear(size_t g, size_t i, const Orientation & Q, const Symmetric &
                    stress) const;

  const std::shared_ptr<SymmetryGroup> symmetry();

 private:
  void make_reciprocal_lattice_();
  static void assert_miller_(std::vector<int> m);

 private:
  Vector a1_, a2_, a3_, b1_, b2_, b3_;
  const std::shared_ptr<SymmetryGroup> symmetry_;
  
  std::vector<std::vector<Vector>> burgers_vectors_;
  std::vector<std::vector<Vector>> slip_directions_;
  std::vector<std::vector<Vector>> slip_planes_;

  std::vector<size_t> offsets_;
};

class CubicLattice: public Lattice {
 public:
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
