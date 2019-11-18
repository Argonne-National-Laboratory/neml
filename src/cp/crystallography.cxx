#include "crystallography.h"

#include "../math/nemlmath.h"

#include <cmath>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include <iostream>

namespace neml {

std::vector<Orientation> symmetry_rotations(std::string sclass)
{
  const double a = sqrt(2.0) / 2.0;
  const double b = sqrt(3.0) / 2.0;
  const double h = 0.5;
  const std::vector<Orientation> tetragonal = {
    Orientation({ 1, 0, 0, 0}),
    Orientation({ 0, 0, 1, 0}),
    Orientation({ 0, 1, 0, 0}),
    Orientation({ 0, 0, 0, 1}),
    Orientation({ a, 0, 0,-a}),
    Orientation({ a, 0, 0, a}),
    Orientation({ 0, a, a, 0}),
    Orientation({ 0,-a, a, 0})
  };
  const std::vector<Orientation> hexagonal = {
    Orientation({ 1, 0, 0, 0}),
    Orientation({ b, 0, 0,-h}),
    Orientation({ b, 0, 0, h}),
    Orientation({-h, 0, 0, b}),
    Orientation({ 0, 0, 0, 1}),
    Orientation({ h, 0, 0, b}),
    Orientation({ 0, b,-h, 0}),
    Orientation({ 0, 1, 0, 0}),
    Orientation({ 0, b, h, 0}),
    Orientation({ 0, h, b, 0}),
    Orientation({ 0, 0, 1, 0}),
    Orientation({ 0,-h, b, 0})
  };
  const std::vector<Orientation> cubic = {
    Orientation({ 1, 0, 0, 0}),
    Orientation({ h, h, h, h}),
    Orientation({-h, h, h, h}),
    Orientation({ h,-h, h, h}),
    Orientation({ h, h,-h, h}),
    Orientation({-h,-h,-h, h}),
    Orientation({ h,-h,-h, h}),
    Orientation({-h,-h, h, h}),
    Orientation({-h, h,-h, h}),
    Orientation({ 0, 0, 1, 0}),
    Orientation({ 0, 0, 0, 1}),
    Orientation({ 0, 1, 0 ,0}),
    Orientation({ 0,-h, 0, h}),
    Orientation({ 0, h, 0, h}),
    Orientation({ h, 0, h, 0}),
    Orientation({ h, 0,-h, 0}),
    Orientation({ 0, 0,-h, h}),
    Orientation({ h, h, 0, 0}),
    Orientation({ h,-h, 0, 0}),
    Orientation({ 0, 0, h, h}),
    Orientation({ 0,-h, h, 0}),
    Orientation({ h, 0, 0,-h}),
    Orientation({ 0, h, h, 0}),
    Orientation({ h, 0, 0, h})
  };

  std::vector<Orientation> res;
  if (sclass == "432") {
    for (int i=0; i<24; i++) {
      res.push_back(cubic[i]);
    }
  }
  else if (sclass == "23") {
    for (int i=0; i<12; i++) {
      res.push_back(cubic[i]);
    }
  }
  else if (sclass == "622") {
    for (int i=0; i<12; i++) {
      res.push_back(hexagonal[i]);
    }
  }
  else if (sclass == "32") {
    for (int i=0; i<3; i++) {
      res.push_back(hexagonal[i]);
    }
    for (int i=9; i<12; i++) {
      res.push_back(hexagonal[i]);
    }
  }
  else if (sclass == "6") {
    for (int i=0; i<6; i++) {
      res.push_back(hexagonal[i]);
    }
  }
  else if (sclass == "3") {
    for (int i=0; i<3; i++) {
      res.push_back(hexagonal[i]);
    }
  }
  else if (sclass == "42") {
    for (int i=0; i<8; i++) {
      res.push_back(tetragonal[i]);
    }
  }
  else if (sclass == "4") {
    res.push_back(tetragonal[0]);
    for (int i=3; i<6; i++) {
      res.push_back(tetragonal[i]);
    }
  }
  else if (sclass == "222") {
    for (int i=0; i<4; i++) {
      res.push_back(tetragonal[i]);
    }
  }
  else if (sclass == "2") {
    for (int i=0; i<2; i++) {
      res.push_back(tetragonal[i]);
    }
  }
  else if (sclass == "1") {
    res.push_back(tetragonal[0]);
  }
  else {
    throw std::invalid_argument("Invalid crystal class " + sclass + ".");
  }

  return res;
}

SymmetryGroup::SymmetryGroup(std::string sclass) :
    ops_(symmetry_rotations(sclass))
{
  misops_.reserve(ops_.size() * ops_.size());
  for (auto a : ops_) {
    for (auto b : ops_) {
      misops_.push_back(a * b);
    }
  }
}

SymmetryGroup::~SymmetryGroup()
{

}

std::string SymmetryGroup::type()
{
  return "SymmetryGroup";
}

ParameterSet SymmetryGroup::parameters()
{
  ParameterSet pset(SymmetryGroup::type());

  pset.add_parameter<std::string>("sclass");

  return pset;
}

std::unique_ptr<NEMLObject> SymmetryGroup::initialize(ParameterSet & params)
{
  return neml::make_unique<SymmetryGroup>(
      params.get_parameter<std::string>("sclass")); 
}

const std::vector<Orientation> & SymmetryGroup::ops() const
{
  return ops_;
}

size_t SymmetryGroup::nops() const
{
  return ops_.size();
}

Orientation SymmetryGroup::misorientation(const Orientation & a, 
                                          const Orientation & b) const
{
  Orientation best;
  Orientation ab = a * b.inverse();
  double angle_best = 2 * M_PI;
  double cn[3];
  double ca;
  for (auto misop : misops_) {
    Orientation trial = misop * ab;
    trial.to_axis_angle(cn, ca);
    if (ca < angle_best) {
      angle_best = ca;
      best = trial;
    }
  }

  return best;
}

std::vector<Orientation> SymmetryGroup::misorientation_block(
    const std::vector<Orientation> & A, const std::vector<Orientation> & B)
{
  size_t N = A.size();
  
  if (N != B.size()) {
    throw std::runtime_error("Blocks of input orientations do not have matching sizes!");
  }

  double * V = new double[N*4];
  for (size_t i = 0; i < N; i++) {
    Orientation v(&V[4*i]);
    v = A[i] * B[i].inverse();
  }
  
  size_t S = misops_.size();
  double * M = new double[S * 4 * 4];
  for (size_t i = 0; i < S; i++) {
    misops_[i].to_product_matrix(&M[i*4*4]);
  }

  // Answer is now M * V.T
  double * R = new double[S * 4 * N];
  mat_mat_ABT(S*4, N, 4, M, V, R);

  delete [] M;
  delete [] V;

  std::vector<Orientation> res(N);
  for (size_t i = 0; i < N; i++) {
    double best_angle = 3 * M_PI;
    size_t bi = -1;
    for (size_t j = 0; j < S; j++) {
      double ang = 2.0 * acos(R[CINDEX((j*4),i,N)]);
      if (ang < best_angle) {
        best_angle = ang;
        bi = j;
      }
    }
    res[i] = Orientation({R[CINDEX((bi*4+0),i,N)],R[CINDEX((bi*4+1),i,N)],R[CINDEX((bi*4+2),i,N)],R[CINDEX((bi*4+3),i,N)]});
  }

  delete [] R;

  return res;
}

Lattice::Lattice(Vector a1, Vector a2, Vector a3,
                 std::shared_ptr<SymmetryGroup> symmetry,
                 list_systems isystems) :
    a1_(a1), a2_(a2), a3_(a3), symmetry_(symmetry), offsets_({0}),
    setup_(false)
{
  make_reciprocal_lattice_();

  for (auto system : isystems) {
    add_slip_system(system.first, system.second);
  }
}

Lattice::~Lattice()
{

}

Vector Lattice::miller2cart_direction(std::vector<int> m)
{
  Lattice::assert_miller_(m);
  std::vector<int> mr = reduce_gcd(m);

  return (double) mr[0] * a1_ + (double) mr[1] * a2_ + (double) mr[2] * a3_;
}

Vector Lattice::miller2cart_plane(std::vector<int> m)
{
  Lattice::assert_miller_(m);
  std::vector<int> mr = reduce_gcd(m);

  return (double) mr[0] * b1_ + (double) mr[1] * b2_ + (double) mr[2] * b3_;
}

std::vector<Vector> Lattice::equivalent_vectors(Vector v)
{
  std::vector<Vector> uniques;
  for (auto R = symmetry_->ops().begin(); R != symmetry_->ops().end(); ++R) {
    Vector curr = R->apply(v);
    bool dup = false;
    for (auto vi = uniques.begin(); vi != uniques.end(); ++vi) {
      if (*vi == curr) {
        dup = true;
        break;
      }
    }
    if (not dup) {
      uniques.push_back(curr);
    }
  }

  return uniques;
}

std::vector<Vector> Lattice::equivalent_vectors_bidirectional(Vector v)
{
  std::vector<Vector> both = equivalent_vectors(v);
  
  std::vector<Vector> single;

  for (auto a = both.begin(); a != both.end(); ++a) {
    bool dup = false;
    for (auto b = single.begin(); b != single.end(); ++b) {
      if (*a == -(*b)) {
        dup = true;
        break;
      }
    }
    if (not dup) {
      single.push_back(*a);
    }
  }
  
  return single;
}

void Lattice::add_slip_system(std::vector<int> d, std::vector<int> p)
{
  std::vector<Vector> burgers, directions, normals;

  std::vector<Vector> pbs = equivalent_vectors_bidirectional(
      miller2cart_direction(d));
  std::vector<Vector> pns = equivalent_vectors_bidirectional(
      miller2cart_plane(p));
  
  for (auto bi = pbs.begin(); bi != pbs.end(); ++bi) {
    for (auto ni = pns.begin(); ni != pns.end(); ++ni) {
      Vector nd = *bi / bi->norm();
      Vector nn = *ni / ni->norm();
      
      if (isclose(nd.dot(nn), 0.0)) {
        burgers.push_back(*bi);
        directions.push_back(nd);
        normals.push_back(nn);
      }
    }
  }

  if (burgers.size() != 0) {
    burgers_vectors_.push_back(burgers);
    slip_directions_.push_back(directions);
    slip_planes_.push_back(normals);
    offsets_.push_back(offsets_.back() + burgers.size());
  }
}

size_t Lattice::ngroup() const
{
  return slip_planes_.size();
}

size_t Lattice::nslip(size_t g) const
{
  return slip_planes_[g].size();
}

size_t Lattice::flat(size_t g, size_t i) const
{
  return offsets_[g] + i;
}

Symmetric Lattice::M(size_t g, size_t i, const Orientation & Q)
{
  cache_rot_(Q);
  return Ms_[g][i];
  return Q.apply(Symmetric(outer(slip_directions_[g][i], 
                                 slip_planes_[g][i])));
}

Skew Lattice::N(size_t g, size_t i, const Orientation & Q)
{
  cache_rot_(Q);
  return Ns_[g][i];
}

double Lattice::shear(size_t g, size_t i, const Orientation & Q, 
                      const Symmetric & stress)
{
  return M(g,i,Q).contract(stress);
}

Symmetric Lattice::d_shear(size_t g, size_t i, const Orientation & Q, 
                           const Symmetric & stress)
{
  return M(g, i, Q);
}

const std::shared_ptr<SymmetryGroup> Lattice::symmetry()
{
  return symmetry_;
}

void Lattice::make_reciprocal_lattice_()
{
  b1_ = a2_.cross(a3_) / a1_.dot(a2_.cross(a3_));
  b2_ = a3_.cross(a1_) / a2_.dot(a3_.cross(a1_));
  b3_ = a1_.cross(a2_) / a3_.dot(a1_.cross(a2_));
}

void Lattice::assert_miller_(std::vector<int> m)
{
  if (m.size() != 3)   {
    throw std::invalid_argument("Miller indices must have length 3");
  }
}

void Lattice::cache_rot_(const Orientation & Q)
{
  if (setup_ and (hash_ == Q.hash())) return;

  setup_ = true;
  hash_ = Q.hash();

  Ms_.resize(ngroup());
  Ns_.resize(ngroup());
  for (size_t g = 0; g < ngroup(); g++) {
    Ms_[g].resize(nslip(g));
    Ns_[g].resize(nslip(g));
    for (size_t i = 0; i < nslip(g); i++) {
      Ms_[g][i] = Q.apply(Symmetric(outer(slip_directions_[g][i], 
                                          slip_planes_[g][i]))); 
      Ns_[g][i] = Q.apply(Skew(outer(slip_directions_[g][i], 
                                     slip_planes_[g][i]))); 
    }
  }
}

CubicLattice::CubicLattice(double a, 
                           list_systems isystems) :
    Lattice(Vector({a,0,0}),Vector({0,a,0}),Vector({0,0,a}), 
            std::make_shared<SymmetryGroup>("432"), isystems)
{

}

std::string CubicLattice::type()
{
  return "CubicLattice";
}

ParameterSet CubicLattice::parameters()
{
  ParameterSet pset(CubicLattice::type());

  pset.add_parameter<double>("a");
  pset.add_optional_parameter<list_systems>("slip_systems", 
                                            list_systems());

  return pset;
}

std::unique_ptr<NEMLObject> CubicLattice::initialize(ParameterSet & params)
{
  return neml::make_unique<CubicLattice>(
      params.get_parameter<double>("a"),
      params.get_parameter<list_systems>("slip_systems")); 
}

} // namespace neml
