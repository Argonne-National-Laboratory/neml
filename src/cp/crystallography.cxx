#include "cp/crystallography.h"

#include "math/nemlmath.h"

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
    Orientation({ -h, 0, 0,b}),
    Orientation({ h, 0, 0, b}),
    Orientation({b, 0, 0, -h}),
    Orientation({ 0, 0, 0, 1}),
    Orientation({ b, 0, 0, h}),
    Orientation({ 0,-h, b, 0}),
    Orientation({ 0, 1, 0, 0}),
    Orientation({ 0, h, b, 0}),
    Orientation({ 0, b, h, 0}),
    Orientation({ 0, 0, 1, 0}),
    Orientation({ 0, b,-h, 0})
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

SymmetryGroup::SymmetryGroup(ParameterSet & params) :
    NEMLObject(params),
    ops_(symmetry_rotations(params.get_parameter<std::string>("sclass")))
{
  misops_.reserve(ops_.size() * ops_.size());
  for (auto a : ops_) {
    for (auto b : ops_) {
      misops_.push_back(a * b);
    }
  }
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
  return neml::make_unique<SymmetryGroup>(params);
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
      double f = std::max(-1.0,std::min(1.0,R[CINDEX((j*4),i,N)])); // need to check to prevent precision error and lie between [-1,1]
      double ang = 2.0 * acos(f);
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

std::shared_ptr<SymmetryGroup> get_group(std::string grp)
{
  ParameterSet params = SymmetryGroup::parameters();
  params.assign_parameter("sclass", grp);

  return std::make_shared<SymmetryGroup>(params);
}

Lattice::Lattice(Vector a1, Vector a2, Vector a3,
                 std::shared_ptr<SymmetryGroup> symmetry,
                 list_systems isystems,
                 twin_systems tsystems) :
    a1_(a1), a2_(a2), a3_(a3), symmetry_(symmetry), 
    offsets_({0}), setup_(false)
{
  make_reciprocal_lattice_();

  for (auto system : isystems) {
    add_slip_system(system.first, system.second);
  }
  for (auto system : tsystems) {
    add_twin_system(std::get<0>(system), std::get<1>(system),
                    std::get<2>(system), std::get<3>(system));
  }
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
  current_slip_.push_back(std::make_pair(d,p));

  std::vector<Vector> burgers, directions, normals;
  std::vector<Orientation> reorientations;

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
        reorientations.push_back(Orientation(std::vector<double>({1,0,0,0})));
      }
    }
  }

  if (burgers.size() != 0) {
    burgers_vectors_.push_back(burgers);
    slip_directions_.push_back(directions);
    slip_planes_.push_back(normals);
    offsets_.push_back(offsets_.back() + burgers.size());
    slip_types_.push_back(Lattice::SlipType::Slip);
    shear_.push_back(0.0);
    reorientations_.push_back(reorientations);
    update_normals_(normals);
  }
}

void Lattice::add_twin_system(std::vector<int> eta1, std::vector<int> K1,
                              std::vector<int> eta2, std::vector<int> K2)
{
  current_twin_.push_back(std::make_tuple(eta1,K1,eta2,K2));

  std::vector<Vector> burgers, directions, normals;
  std::vector<Orientation> reorientations;

  std::vector<Vector> cart_K1 = equivalent_vectors(miller2cart_plane(K1));
  std::vector<Vector> cart_eta2 = equivalent_vectors(miller2cart_direction(eta2));

  std::vector<Vector> cart_eta1 =
      equivalent_vectors(miller2cart_direction(eta1));
  std::vector<Vector> cart_K2 = equivalent_vectors(miller2cart_plane(K2));
  
  // Keep this available for the group
  double shear = 0;

  for (auto K1i = cart_K1.begin(); K1i != cart_K1.end(); ++K1i) {
    for (auto eta1i = cart_eta1.begin(); eta1i != cart_eta1.end(); ++eta1i) {
      for (auto K2i = cart_K2.begin(); K2i != cart_K2.end(); ++K2i) {
        for (auto eta2i = cart_eta2.begin(); eta2i != cart_eta2.end(); ++eta2i) {
          // Make unit vectors
          Vector l = *eta1i / eta1i->norm();
          Vector m = *K1i / K1i->norm();
          Vector n = *K2i / K2i->norm();
          Vector g = *eta2i / eta2i->norm();
          
          // Check geometry of the shearing direction from K1 and eta2
          double val = g.dot(m);
          if (isclose(val,0)) continue;
          Vector eta1_tr = 2*(m-g/val);
          Vector l_tr = eta1_tr/eta1_tr.norm();
          if (! isclose(fabs(l_tr.dot(l)), 1)) continue;

          // Check the geometry of the undisturbed plane (use g and m and l...)
          Vector P = l.cross(m);
          P /= P.norm();
          Vector n_tr = g.cross(P);
          n_tr /= n_tr.norm();
          if (! isclose(fabs(n_tr.dot(n)), 1)) continue;

          // Requirements:
          // 3) Angle between l and g is obtuse
          // 4) Angle between l and n is acute
          // 5) Angle between g and m is acute
          // The 0.5 is some weird thing going on with Ti HCP Compression Twin
          if ((l.dot(g) < 0) && (l.dot(n) > 0.5) &&
              (g.dot(m) > 0.5)) {

            // Check for the degenerate case
            bool duplicate = false;
            for (size_t i = 0; i < directions.size(); i++) {
              if ((isclose(directions[i].dot(l), -1) &&
                  isclose(normals[i].dot(m), -1)) ||
                  (isclose(directions[i].dot(l),1) &&
                   isclose(normals[i].dot(m),1))) {
                duplicate = true;
                break;
              }
            }
            if (duplicate) continue;

            // Right for compound twin
            auto q = Orientation::createAxisAngle(m.data(), 180, "degrees");
            reorientations.push_back(q);

            // Shear
            shear = sqrt(4.0*(1.0/std::pow(g.dot(m),2)-1.0));

            burgers.push_back(*eta1i);
            directions.push_back(l);
            normals.push_back(m);
          }
        }
      }
    }
  }
  
  if (burgers.size() != 0) {
    burgers_vectors_.push_back(burgers);
    slip_directions_.push_back(directions);
    slip_planes_.push_back(normals);
    offsets_.push_back(offsets_.back() + burgers.size());
    slip_types_.push_back(Lattice::SlipType::Twin);
    shear_.push_back(shear);
    reorientations_.push_back(reorientations);
    update_normals_(normals);
  }
}

size_t Lattice::ntotal() const
{
  size_t total = 0;
  
  for (size_t i = 0; i < ngroup(); i++) {
    total += nslip(i);
  }

  return total;
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

Lattice::SlipType Lattice::slip_type(size_t g, size_t i) const
{
  return slip_types_[g];
}

double Lattice::characteristic_shear(size_t g, size_t i) const
{
  return shear_[g];
}

Orientation Lattice::reorientation(size_t g, size_t i) const
{
  return reorientations_[g][i];
}

double Lattice::burgers(size_t g, size_t i) const
{
  return burgers_vectors_[g][i].norm();
}

const Symmetric & Lattice::M(size_t g, size_t i, const Orientation & Q)
{
  cache_rot_(Q);
  return Ms_[g][i];
}

const Skew & Lattice::N(size_t g, size_t i, const Orientation & Q)
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

const std::vector<Vector> Lattice::unique_planes() const
{
  return normals_;
}

size_t Lattice::nplanes() const
{
  return normals_.size();
}

size_t Lattice::plane_index(size_t g, size_t i) const
{
  return normal_map_[g][i];
}

std::vector<std::pair<size_t,size_t>> Lattice::plane_systems(size_t i) const
{
  std::vector<std::pair<size_t,size_t>> res;
  for (size_t g = 0; g < ngroup(); g++) {
    for (size_t j = 0; j < nslip(g); j++) {
      if (plane_index(g,j) == i) res.push_back(std::make_pair(g,j));
    }
  }
  return res;
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

void Lattice::update_normals_(const std::vector<Vector> & new_planes)
{
  std::vector<size_t> new_indices;

  for (size_t i = 0; i < new_planes.size(); i++) {
    size_t j = 0;
    for (; j < normals_.size(); j++) {
      if (isclose(fabs(normals_[j].dot(new_planes[i])), 1)) {
        new_indices.push_back(j);
        break;
      }
    }
    if (j == normals_.size()) {
      normals_.push_back(new_planes[i]);
      new_indices.push_back(normals_.size()-1);
    }
  }
  normal_map_.push_back(new_indices);
}

GeneralLattice::GeneralLattice(ParameterSet & params) :
    NEMLObject(params),
    Lattice(
        Vector(params.get_parameter<std::vector<double>>("a1")),
        Vector(params.get_parameter<std::vector<double>>("a2")),
        Vector(params.get_parameter<std::vector<double>>("a3")),
        params.get_object_parameter<SymmetryGroup>("symmetry_group"),
        params.get_parameter<list_systems>("slip_systems"), 
        params.get_parameter<twin_systems>("twin_systems"))
{

}

std::string GeneralLattice::type()
{
  return "GeneralLattice";
}

ParameterSet GeneralLattice::parameters()
{
  ParameterSet pset(GeneralLattice::type());

  pset.add_parameter<std::vector<double>>("a1");
  pset.add_parameter<std::vector<double>>("a2");
  pset.add_parameter<std::vector<double>>("a3");
  pset.add_parameter<NEMLObject>("symmetry_group");
  pset.add_optional_parameter<list_systems>("slip_systems",
                                            list_systems());
  pset.add_optional_parameter<list_systems>("twin_systems",
                                            twin_systems());

  return pset;
}

ParameterSet & GeneralLattice::current_parameters()
{
  current_params_.assign_parameter("slip_systems", current_slip_);
  current_params_.assign_parameter("twin_systems", current_twin_);
  return current_params_;
}

std::unique_ptr<NEMLObject> GeneralLattice::initialize(ParameterSet & params)
{
  return neml::make_unique<GeneralLattice>(params);
}

CubicLattice::CubicLattice(ParameterSet & params) :
    NEMLObject(params),
    Lattice(
        Vector({params.get_parameter<double>("a"),0,0}),
        Vector({0,params.get_parameter<double>("a"),0}),
        Vector({0,0,params.get_parameter<double>("a")}),
        get_group("432"), 
        params.get_parameter<list_systems>("slip_systems"), 
        params.get_parameter<twin_systems>("twin_systems"))
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
  pset.add_optional_parameter<list_systems>("twin_systems",
                                            twin_systems());

  return pset;
}

ParameterSet & CubicLattice::current_parameters()
{
  current_params_.assign_parameter("slip_systems", current_slip_);
  current_params_.assign_parameter("twin_systems", current_twin_);
  return current_params_;
}

std::unique_ptr<NEMLObject> CubicLattice::initialize(ParameterSet & params)
{
  return neml::make_unique<CubicLattice>(params);
}

HCPLattice::HCPLattice(ParameterSet & params) :
    NEMLObject(params),
    Lattice(
        Vector({params.get_parameter<double>("a")/2,-sqrt(3)*params.get_parameter<double>("a")/2,0}),
        Vector({params.get_parameter<double>("a")/2,sqrt(3)*params.get_parameter<double>("a")/2,0}),
        Vector({0,0,params.get_parameter<double>("c")}),
        get_group("622"))
{
  list_systems isystems = params.get_parameter<list_systems>("slip_systems"); 
  twin_systems tsystems = params.get_parameter<twin_systems>("twin_systems");


  for (auto system : isystems) {
    add_slip_system(system.first, system.second);
  }
  for (auto system : tsystems) {
    add_twin_system(std::get<0>(system), std::get<1>(system),
                    std::get<2>(system), std::get<3>(system));
  }
}

std::string HCPLattice::type()
{
  return "HCPLattice";
}

ParameterSet HCPLattice::parameters()
{
  ParameterSet pset(HCPLattice::type());

  pset.add_parameter<double>("a");
  pset.add_parameter<double>("c");
  pset.add_optional_parameter<list_systems>("slip_systems",
                                            list_systems());
  pset.add_optional_parameter<twin_systems>("twin_systems",
                                            twin_systems());

  return pset;
}

std::unique_ptr<NEMLObject> HCPLattice::initialize(ParameterSet & params)
{
  return neml::make_unique<HCPLattice>(params);
}

Vector HCPLattice::miller2cart_direction(std::vector<int> m)
{
  assert_miller_bravais_(m);

  return Lattice::miller2cart_direction({2*m[0] + m[1], 2*m[1] + m[0], m[3]});
}

Vector HCPLattice::miller2cart_plane(std::vector<int> m)
{
  assert_miller_bravais_(m);
  return Lattice::miller2cart_plane({m[0], m[1], m[3]});
}

ParameterSet & HCPLattice::current_parameters()
{
  current_params_.assign_parameter("slip_systems", current_slip_);
  current_params_.assign_parameter("twin_systems", current_twin_);
  return current_params_;
}

void HCPLattice::assert_miller_bravais_(std::vector<int> m)
{
  if (m.size() != 4)   
    throw std::invalid_argument("Miller-Bravais indices must have length 4");
  

  if (m[0] + m[1] + m[2] != 0)
    throw std::invalid_argument("h + k + i must sum to 0 for Miller-Bravais"
                                " notation");
}

} // namespace neml
