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

NEML_EXPORT std::vector<Orientation> symmetry_rotations(std::string sclass);

class NEML_EXPORT SymmetryGroup: public NEMLObject {
 public:
  /// Initialize with the Hermann-Mauguin notation as a string
  SymmetryGroup(ParameterSet & params);

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

  /// Find the disorientation in a blocked way that trades memory for cpu
  std::vector<Orientation> misorientation_block(const std::vector<Orientation> & A, const std::vector<Orientation> & B);

 private:
  const std::vector<Orientation> ops_;
  std::vector<Orientation> misops_;
};

static Register<SymmetryGroup> regSymmetryGroup;

std::shared_ptr<SymmetryGroup> get_group(std::string);

class NEML_EXPORT Lattice {
 public:
  /// Initialize with the three lattice vectors, the symmetry group and
  /// (optionally) a initial list of slip systems
  Lattice(Vector a1, Vector a2, Vector a3, std::shared_ptr<SymmetryGroup> symmetry,
          list_systems isystems = {}, twin_systems tsystems = {});
  virtual ~Lattice() {}; // clang??

  /// Type: slip or twin
  enum SlipType {Slip=0, Twin=1};

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
  /// Return the list of normalized slip normals
  const std::vector<std::vector<Vector>> & slip_planes() {return slip_planes_;};
  /// Return the list of characteristic shears
  const std::vector<double> characteristic_shears() {return shear_;}; 
  /// Return the list of slip system types
  const std::vector<SlipType> & slip_types() {return slip_types_;};

  /// Convert Miller directions to cartesian vectors
  virtual Vector miller2cart_direction(std::vector<int> m);
  /// Convert Miller planes to cartesian normal vectors
  virtual Vector miller2cart_plane(std::vector<int> m);

  /// Find all sets of equivalent vectors (+/- different)
  std::vector<Vector> equivalent_vectors(Vector v);
  /// Find all all sets of equivalent vectors (+/- the same)
  std::vector<Vector> equivalent_vectors_bidirectional(Vector v);

  /// Add a slip system given the Miller direction and plane
  virtual void add_slip_system(std::vector<int> d, std::vector<int> p);

  /// Add a twin system given the twin direction and plane
  virtual void add_twin_system(std::vector<int> eta1, std::vector<int> K1,
                               std::vector<int> eta2, std::vector<int> K2);

  /// Number of total slip systems
  size_t ntotal() const;
  /// Number of groups of slip systems
  size_t ngroup() const;
  /// Number of slip systems in group g
  size_t nslip(size_t g) const;
  /// Flat index of slip group g, system i
  size_t flat(size_t g, size_t i) const;
  
  /// Type of system: slip or twin
  SlipType slip_type(size_t g, size_t i) const;
  /// Characteristic shear (0 for slip systems)
  double characteristic_shear(size_t g, size_t i) const;
  /// Twin reorientation operator (identity for slip systems)
  Orientation reorientation(size_t g, size_t i) const;

  /// Norm of the Burgers vector for a particular system
  double burgers(size_t g, size_t i) const;

  /// Return the sym(d x n) tensor for group g, system i, rotated with Q
  const Symmetric & M(size_t g, size_t i, const Orientation & Q);
  /// Return the skew(d x n) tensor for group g, system i, rotated with Q
  const Skew & N(size_t g, size_t i, const Orientation & Q);

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

  /// Return a list of Cartesian vectors giving the unique slip planes
  const std::vector<Vector> unique_planes() const;
  /// The number of unique slip planes
  size_t nplanes() const;
  /// Given a slip system return the index into the unique slip planes
  size_t plane_index(size_t g, size_t i) const;
  /// Given the unique slip plane index return the vector of (g,i) tuples
  std::vector<std::pair<size_t,size_t>> plane_systems(size_t i) const;

 private:
  void make_reciprocal_lattice_();
  static void assert_miller_(std::vector<int> m);

  void cache_rot_(const Orientation & Q);

  void update_normals_(const std::vector<Vector> & new_planes);

 protected:
  list_systems current_slip_;
  twin_systems current_twin_;

 private:
  Vector a1_, a2_, a3_, b1_, b2_, b3_;
  const std::shared_ptr<SymmetryGroup> symmetry_;

  std::vector<std::vector<Vector>> burgers_vectors_;
  std::vector<std::vector<Vector>> slip_directions_;
  std::vector<std::vector<Vector>> slip_planes_;
  std::vector<SlipType> slip_types_;
  std::vector<double> shear_;
  std::vector<std::vector<Orientation>> reorientations_;

  std::vector<size_t> offsets_;

  // Used for caching common asks
  bool setup_;
  size_t hash_;
  std::vector<std::vector<Symmetric>> Ms_;
  std::vector<std::vector<Skew>> Ns_;

  // Used for the normal damage system
  std::vector<Vector> normals_;
  std::vector<std::vector<size_t>> normal_map_;
};

/// General lattice with no specialization
class NEML_EXPORT GeneralLattice: public NEMLObject, public Lattice {
 public:
  GeneralLattice(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Override serialization to account for dynamic changes
  virtual ParameterSet & current_parameters();
};

static Register<GeneralLattice> regGeneralLattice;

class NEML_EXPORT CubicLattice: public NEMLObject, public Lattice {
 public:
  /// Specialized Lattice for cubic systems, initialize with the lattice
  /// parameter
  CubicLattice(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Override serialization to account for dynamic changes
  virtual ParameterSet & current_parameters();
};

static Register<CubicLattice> regCubicLattice;

class NEML_EXPORT HCPLattice: public NEMLObject, public Lattice {
 public:
  /// Specialized to HCP, initialize with a and c
  HCPLattice(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Default parameters
  static ParameterSet parameters();

  /// Convert Miller-Bravais directions to cartesian vectors
  virtual Vector miller2cart_direction(std::vector<int> m);
  /// Convert Miller-Bravais planes to cartesian normal vectors
  virtual Vector miller2cart_plane(std::vector<int> m);

  /// Override serialization to account for dynamic changes
  virtual ParameterSet & current_parameters();

 private:
  void assert_miller_bravais_(std::vector<int> m);
};

static Register<HCPLattice> regHCPLattice;

} // namespace neml

#endif // CRYSTALLOGRAPHY_H
