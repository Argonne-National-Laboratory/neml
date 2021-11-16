#include "math/rotations.h"

#include "math/nemlmath.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <chrono>
#include <functional>

namespace neml {

Quaternion::Quaternion()
{
  alloc_();
  quat_[0] = 0.0;
  quat_[1] = 0.0;
  quat_[2] = 0.0;
  quat_[3] = 1.0;
}

Quaternion::Quaternion(const std::vector<double> v)
{
  alloc_();
  std::copy(v.begin(), v.end(), quat_);
}

Quaternion::Quaternion(double * v)
{
  store_ = false;
  quat_ = v;
}

Quaternion::Quaternion(const Quaternion & other) :
    store_(other.store())
{
  if (other.store()) {
    alloc_();
    std::copy(other.quat(), other.quat()+4, quat_);
  }
  else {
    quat_ = const_cast<double*>(other.quat());
  }
}

Quaternion::Quaternion(const Quaternion && other) :
    store_(other.store())
{
  if (store_) {
    alloc_();
    std::copy(other.quat(), other.quat()+4, quat_);
  }
  else {
    quat_ = const_cast<double*>(other.quat());
  }
}

Quaternion::~Quaternion()
{
  if (store_) {
    delete [] quat_;
  }
  quat_ = nullptr;
}

Quaternion & Quaternion::operator=(const Quaternion & rhs)
{
  // Copy
  if (this != &rhs) {
    std::copy(rhs.quat(), rhs.quat()+4, quat_);
  };
  return *this;
}

Quaternion & Quaternion::operator=(const Quaternion && rhs)
{
  std::copy(rhs.quat(), rhs.quat()+4, quat_);

  return *this;
}

bool Quaternion::store() const
{
  return store_;
}

const double * Quaternion::quat() const
{
  return quat_;
}

double * Quaternion::data()
{
  return quat_;
}

double Quaternion::norm() const
{
  return sqrt(
      quat_[0] * quat_[0]
    + quat_[1] * quat_[1]
    + quat_[2] * quat_[2]
    + quat_[3] * quat_[3]
    );
}

Quaternion Quaternion::opposite() const
{
  Quaternion nq;
  opposite_(nq.data());

  return nq;
}

Quaternion Quaternion::operator-() const
{
  return this->opposite();
}

Quaternion Quaternion::conj() const
{
  Quaternion nq;
  conj_(nq.data());

  return nq;
}

Quaternion Quaternion::flip() const
{
  Quaternion nq;
  flip_(nq.data());

  return nq;
}

Quaternion Quaternion::inverse() const
{
  Quaternion nq;
  inverse_(nq.data());

  return nq;
}

Quaternion Quaternion::exp() const
{
  Quaternion nq;
  exp_(nq.data());

  return nq;
}

Quaternion Quaternion::log() const
{
  Quaternion nq;
  log_(nq.data());

  return nq;
}

Quaternion & Quaternion::operator*=(const Quaternion & rhs)
{
  multiply_(rhs.quat());

  return *this;
}

Quaternion & Quaternion::operator*=(double scalar)
{
  smultiply_(scalar);

  return *this;
}

Quaternion & Quaternion::operator/=(const Quaternion & rhs)
{
  Quaternion inv = rhs.inverse();
  return this->operator*=(inv);
}

Quaternion & Quaternion::operator/=(double scalar)
{
  return this->operator*=(1.0 / scalar);
}

Quaternion Quaternion::pow(double w) const
{
  return (w * this->log()).exp();
}

double Quaternion::dot(const Quaternion & other) const
{
  double s = 0.0;
  for (int i=0; i<4; i++) {
    s += quat()[i] * other.quat()[i];
  }
  return s;
}

size_t Quaternion::hash() const
{
  size_t key = 0;

  for (size_t i = 0; i<4; i++) {
    // One-liner from boost
    key ^= (std::hash<double>{}(quat_[i]) + 0x9e3779b9 + (key<<6) + (key>>2));
  }

  return key;
}

void Quaternion::to_product_matrix(double * M) const
{
  M[0] = quat_[0];
  M[1] = -quat_[1];
  M[2] = -quat_[2];
  M[3] = -quat_[3];

  M[4] = quat_[1];
  M[5] = quat_[0];
  M[6] = -quat_[3];
  M[7] = quat_[2];

  M[8] = quat_[2];
  M[9] = quat_[3];
  M[10] = quat_[0];
  M[11] = -quat_[1];

  M[12] = quat_[3];
  M[13] = -quat_[2];
  M[14] = quat_[1];
  M[15] = quat_[0];
}

Quaternion operator*(double s, const Quaternion & q)
{
  Quaternion cpy(q);
  cpy *= s;
  return cpy;
}

Quaternion operator*(const Quaternion & q, double s)
{
  return s * q;
}

Quaternion operator*(const Quaternion & lhs, const Quaternion & rhs)
{
  Quaternion cpy(lhs);
  cpy *= rhs;
  return cpy;
}

Quaternion operator/(const Quaternion & q, double s)
{
  Quaternion cpy(q);
  cpy /= s;
  return cpy;
}

Quaternion operator/(const Quaternion & lhs, const Quaternion & rhs)
{
  Quaternion cpy(lhs);
  cpy /= rhs;
  return cpy;
}

// Actual math definitions
void Quaternion::opposite_(double * const out) const
{
  out[0] = -quat_[0];
  out[1] = -quat_[1];
  out[2] = -quat_[2];
  out[3] = -quat_[3];
}

void Quaternion::conj_(double * const out) const
{
  out[0] = quat_[0];
  out[1] = -quat_[1];
  out[2] = -quat_[2];
  out[3] = -quat_[3];
}

void Quaternion::flip_(double * const out) const
{
  out[0] = -quat_[0];
  out[1] = quat_[1];
  out[2] = quat_[2];
  out[3] = quat_[3];
}

void Quaternion::inverse_(double * const out) const
{
  double ns = this->norm();
  ns *= ns;
  out[0] = quat_[0] / ns;
  out[1] = -quat_[1] / ns;
  out[2] = -quat_[2] / ns;
  out[3] = -quat_[3] / ns;
}

void Quaternion::exp_(double * const out) const
{
  double nv = 0.0;
  for (int i=1; i<4; i++)
  {
    nv += quat_[i] * quat_[i];
  }
  nv = sqrt(nv);

  out[0] = cos(nv);
  double sf = sin(nv) / nv;
  out[1] = quat_[1] * sf;
  out[2] = quat_[2] * sf;
  out[3] = quat_[3] * sf;
}

void Quaternion::log_(double * const out) const
{
  double n = this->norm();
  double nv = 0.0;
  for (int i=1; i<4; i++)
  {
    nv += quat_[i] * quat_[i];
  }
  nv = sqrt(nv);

  out[0] = ::log(n);
  double sf = acos(quat_[0] / n);
  out[1] = quat_[1] / nv * sf;
  out[2] = quat_[2] / nv * sf;
  out[3] = quat_[3] / nv * sf;
}

void Quaternion::multiply_(const double * const b)
{
  double c[4];
  c[0] = quat_[0]*b[0] - (quat_[1]*b[1] + quat_[2]*b[2] + quat_[3]*b[3]);
  c[1] = quat_[0]*b[1] + quat_[1]*b[0] + quat_[2]*b[3] - quat_[3]*b[2];
  c[2] = quat_[0]*b[2] + quat_[2]*b[0] + quat_[3]*b[1] - quat_[1]*b[3];
  c[3] = quat_[0]*b[3] + quat_[3]*b[0] + quat_[1]*b[2] - quat_[2]*b[1];
  std::copy(c, c+4, quat_);
}

void Quaternion::smultiply_(double s)
{
  for (int i=0; i<4; i++) {
    quat_[i] *= s;
  }
}

void Quaternion::alloc_()
{
  store_ = true;
  quat_ = new double[4];
}

std::ostream & operator<<(std::ostream & os, const Quaternion & q)
{
  const double * const qv = q.quat();
  os << "[";
  for (size_t i = 0; i < 4; i++) {
    os << qv[i] << " ";
  }
  os << "]";

  return os;
}

// Unit quaternion stuff
Orientation Orientation::createRodrigues(const double * const r)
{
  Orientation q;
  q.setRodrigues(r);
  return q;
}

void Orientation::setRodrigues(const double * const r)
{
  double nr = 0.0;
  for (int i=0; i<3; i++) {
    nr += r[i] * r[i];
  }

  double f = 1.0 / (sqrt(1.0 + nr));

  quat_[0] = f;
  quat_[1] = r[0] * f;
  quat_[2] = r[1] * f;
  quat_[3] = r[2] * f;
}

Orientation Orientation::createMatrix(const double * const M)
{
  Orientation q;
  q.setMatrix(M);
  return q;
}

void Orientation::setMatrix(const double * const M)
{
  double tr = 0.0;
  for (int i=0; i<3; i++) {
    tr += M[CINDEX(i,i, 3)];
  }

  double S, s, x, y, z;
  if (tr > 0.0) {
    S = sqrt(1.0 + tr) * 2.0;
    s = 0.25 * S;
    x = (M[CINDEX(2,1,3)] - M[CINDEX(1,2,3)])/S;
    y = (M[CINDEX(0,2,3)] - M[CINDEX(2,0,3)])/S;
    z = (M[CINDEX(1,0,3)] - M[CINDEX(0,1,3)])/S;
  }
  else if ((M[CINDEX(0,0,3)] > M[CINDEX(1,1,3)]) && (M[CINDEX(0,0,3)] > M[CINDEX(2,2,3)])) {
    S = sqrt(1.0 + M[CINDEX(0,0,3)] - M[CINDEX(1,1,3)] - M[CINDEX(2,2,3)]) * 2.0;
    s = (M[CINDEX(2,1,3)] - M[CINDEX(1,2,3)]) / S;
    x = 0.25*S;
    y = (M[CINDEX(0,1,3)] + M[CINDEX(1,0,3)]) / S;
    z = (M[CINDEX(0,2,3)] + M[CINDEX(2,0,3)])/ S;
  }
  else if (M[CINDEX(1,1,3)] > M[CINDEX(2,2,3)]) {
    S = sqrt(1.0 + M[CINDEX(1,1,3)] - M[CINDEX(0,0,3)] - M[CINDEX(2,2,3)]) * 2.0;
    s = (M[CINDEX(0,2,3)] - M[CINDEX(2,0,3)]) / S;
    x = (M[CINDEX(0,1,3)] + M[CINDEX(1,0,3)]) / S;
    y = 0.25 * S;
    z = (M[CINDEX(1,2,3)] + M[CINDEX(2,1,3)]) / S;
  }
  else {
    S = sqrt(1.0 + M[CINDEX(2,2,3)] - M[CINDEX(0,0,3)] - M[CINDEX(1,1,3)]) * 2.0;
    s = (M[CINDEX(1,0,3)] - M[CINDEX(0,1,3)]) / S;
    x = (M[CINDEX(0,2,3)] + M[CINDEX(2,0,3)]) / S;
    y = (M[CINDEX(1,2,3)] + M[CINDEX(2,1,3)]) / S;
    z = 0.25 * S;
  }

  quat_[0] = s;
  quat_[1] = x;
  quat_[2] = y;
  quat_[3] = z;
}

Orientation Orientation::createAxisAngle(const double * const n, double a,
                                         std::string angles)
{
  Orientation q;
  q.setAxisAngle(n, a, angles);
  return q;
}

void Orientation::setAxisAngle(const double * const n, double a,
                               std::string angles)
{
  a = convert_angle(a, angles);

  quat_[0] = cos(a/2.0);
  double s = sin(a/2.0);
  quat_[1] = s * n[0];
  quat_[2] = s * n[1];
  quat_[3] = s * n[2];
}

Orientation Orientation::createEulerAngles(double a, double b, double c,
                                           std::string angles,
                                           std::string convention)
{
  Orientation q;
  q.setEulerAngles(a, b, c, angles, convention);
  return q;
}

void Orientation::setEulerAngles(double a, double b, double c,
                                 std::string angles,
                                 std::string convention)
{
  a = convert_angle(a, angles);
  b = convert_angle(b, angles);
  c = convert_angle(c, angles);
  double na, nb, nc;
  to_kocks_(a, b, c, na, nb, nc, convention);

  double M[9];
  kocks_to_matrix_(na, nb, nc, M);

  setMatrix(M);
}

Orientation Orientation::createHopf(double alpha, double beta, double gamma,
                                    std::string angles)
{
  Orientation q;
  q.setHopf(alpha, beta, gamma, angles);
  return q;
}

// alpha in [0,2pi], beta in [0,pi] (sphere), gamma in [0,2pi] (circle)
void Orientation::setHopf(double alpha, double beta, double gamma,
                          std::string angles)
{
  alpha = convert_angle(alpha, angles);
  beta = convert_angle(beta, angles);
  gamma = convert_angle(gamma, angles);

  quat_[0] = cos(beta / 2.0) * cos(gamma / 2.0);
  quat_[1] = cos(beta / 2.0) * sin(gamma / 2.0);
  quat_[2] = sin(beta / 2.0) * cos(alpha + gamma / 2.0);
  quat_[3] = sin(beta / 2.0) * sin(alpha + gamma / 2.0);
}

Orientation Orientation::createHyperspherical(double a1, double a2, double a3,
                                              std::string angles)
{
  Orientation q;
  q.setHyperspherical(a1, a2, a3, angles);
  return q;
}

void Orientation::setHyperspherical(double a1, double a2, double a3,
                                    std::string angles)
{
  a1 = convert_angle(a1, angles);
  a2 = convert_angle(a2, angles);
  a3 = convert_angle(a3, angles);

  quat_[0] = cos(a1);
  quat_[1] = sin(a1)*cos(a2);
  quat_[2] = sin(a1)*sin(a2)*cos(a3);
  quat_[3] = sin(a1)*sin(a2)*sin(a3);
}

Orientation Orientation::createVectors(const Vector & x, const Vector & y)
{
  Orientation q;
  q.setVectors(x,y);
  return q;
}

void Orientation::setVectors(const Vector & x, const Vector & y)
{
  if (fabs(x.dot(y)) > 1.0e-16) {
    throw std::runtime_error("Input vectors are not orthonormal!");
  }
  Vector z = x.cross(y);
  double M[9];

  M[0] = x(0);
  M[1] = y(0);
  M[2] = z(0);
  M[3] = x(1);
  M[4] = y(1);
  M[5] = z(1);
  M[6] = x(2);
  M[7] = y(2);
  M[8] = z(2);

  setMatrix(M);
}

Orientation::Orientation() :
    Quaternion()
{

}

Orientation::Orientation(double * v) :
    Quaternion(v)
{
  normalize_();
}

Orientation::Orientation(const std::vector<double> v) :
    Quaternion(v)
{
  normalize_();
}

Orientation::Orientation(const Quaternion & other) :
    Quaternion(other)
{
  normalize_();
}

Orientation Orientation::deepcopy() const
{
  Orientation other;
  
  std::copy(quat_, quat_+4, other.data());

  return other;
}

void Orientation::to_euler(double & a, double & b, double & c,
                           std::string angles, std::string convention) const
{
  double tol = 1.0e-15;
  double M[9];
  to_matrix(M);
  double sa, sb, sc;
  if (fabs(fabs(M[CINDEX(2,2,3)]) - 1.0) < tol) {
    sb = 0.0;
    sa = atan2(M[CINDEX(0,1,3)],M[CINDEX(0,0,3)]) / 2.0 - M_PI/2.0;
    sc = -sa;
  }
  else {
    sb = acos(M[CINDEX(2,2,3)]);
    double s = sin(sb);
    sc = M_PI / 2 - atan2(M[CINDEX(0,2,3)] / s, M[CINDEX(1,2,3)] / s);
    sa = M_PI / 2 - atan2(M[CINDEX(2,0,3)] / s, M[CINDEX(2,1,3)] / s);
  }

  double na, nb, nc;
  from_kocks_(sa, sb, sc, na, nb, nc, convention);
  a = cast_angle(na, angles);
  b = cast_angle(nb, angles);
  c = cast_angle(nc, angles);
}

void Orientation::to_axis_angle(double * const n, double & a,
                                std::string angles) const
{
  double f = std::max(-1.0,std::min(quat_[0],1.0)); // need to check to prevent precision error of > 1
  double na = 2.0 * std::acos(f);
  a = cast_angle(na, angles);

  double s = sin(na / 2.0);
  if (a < 1.0e-16) {  // obvious hack
    // If the angle is zero the axis is arbitrary
    n[0] = 1.0;
    n[1] = 0.0;
    n[2] = 0.0;
  }
  else {
    n[0] = quat_[1] / s;
    n[1] = quat_[2] / s;
    n[2] = quat_[3] / s;
  }
}

void Orientation::to_matrix(double * const M) const
{
  double v1s = quat_[1] * quat_[1];
  double v2s = quat_[2] * quat_[2];
  double v3s = quat_[3] * quat_[3];

  M[0] = 1-2*v2s - 2*v3s;
  M[1] = 2*(quat_[1]*quat_[2] - quat_[3]*quat_[0]);
  M[2] = 2*(quat_[1]*quat_[3] + quat_[2]*quat_[0]);
  M[3] = 2*(quat_[1]*quat_[2] + quat_[3]*quat_[0]);
  M[4] = 1-2*v1s - 2*v3s;
  M[5] = 2*(quat_[2]*quat_[3] - quat_[1]*quat_[0]);
  M[6] = 2*(quat_[1]*quat_[3] - quat_[2]*quat_[0]);
  M[7] = 2*(quat_[2]*quat_[3] + quat_[1]*quat_[0]);
  M[8] = 1-2*v1s-2*v2s;
}

RankTwo Orientation::to_tensor() const
{
  RankTwo res;
  to_matrix(res.s());
  return res;
}

void Orientation::to_rodrigues(double * const v) const
{
  double f;
  if (quat_[0] == 0.0) {
    f = 1.0e-16; // Obvious problem
  }
  else {
    f = quat_[0];
  }

  v[0] = quat_[1] / f;
  v[1] = quat_[2] / f;
  v[2] = quat_[3] / f;
}

// alpha = [0,2pi], beta = [0,pi], gamma = [0,2pi]
// sphere then circle
void Orientation::to_hopf(double & alpha, double & beta, double & gamma,
                          std::string angles) const
{
  double x1 = quat_[0];
  double x2 = quat_[1];
  double x3 = quat_[2];
  double x4 = quat_[3];

  gamma = 2.0 * atan2(x2, x1);
  alpha = atan2(x4, x3) - gamma / 2.0;
  beta = 2.0 * asin(sqrt(x3*x3+x4*x4));

  alpha = cast_angle(alpha, angles);
  beta = cast_angle(beta, angles);
  gamma = cast_angle(gamma, angles);
}

// Last angle has the [0,2pi] range
void Orientation::to_hyperspherical(double & a1, double & a2, double & a3,
                                    std::string angles) const
{
  double x1 = quat_[0];
  double x2 = quat_[1];
  double x3 = quat_[2];
  double x4 = quat_[3];

  a1 = acos(x1 / sqrt(x4*x4 + x3*x3 + x2*x2 + x1*x1));
  a2 = acos(x2 / sqrt(x4*x4 + x3*x3 + x2*x2));
  if (x4 >= 0.0) {
    a3 = acos(x3 / sqrt(x4*x4 + x3*x3));
  }
  else {
    a3 = -acos(x3 / sqrt(x4*x4 + x3*x3));
  }

  a1 = cast_angle(a1, angles);
  a2 = cast_angle(a2, angles);
  a3 = cast_angle(a3, angles);
}

Orientation Orientation::opposite() const
{
  Orientation nq;
  opposite_(nq.data());

  return nq;
}

Orientation Orientation::operator-() const
{
  return this->opposite();
}

Orientation Orientation::conj() const
{
  Orientation nq;
  conj_(nq.data());

  return nq;
}

Orientation Orientation::flip() const
{
  Orientation nq;
  flip_(nq.data());

  return nq;
}

Orientation Orientation::inverse() const
{
  Orientation nq;
  conj_(nq.data());

  return nq;
}

Orientation & Orientation::operator*=(const Orientation & rhs)
{
  multiply_(rhs.quat());

  return *this;
}

Orientation & Orientation::operator/=(const Orientation & rhs)
{
  Orientation inv = rhs.inverse();
  return this->operator*=(inv);
}

Orientation Orientation::pow(double w) const
{
  return Orientation(((w * this->log()).exp()));
}

void Orientation::normalize_()
{
  double nv = this->norm();
  // I'm debating not handling this case...
  if (nv == 0) {
    for (int i=0; i<3; i++) {
      quat_[i] = 0.0;
    }
    quat_[3] = 1.0;
    return;
  }
  for (int i=0; i<4; i++) {
    quat_[i] /= nv;
  }
}

void Orientation::to_kocks_(double a, double b, double c, double & oa,
                            double & ob, double & oc, std::string type)
{
  if (std::string("kocks").compare(type) == 0) {
    oa = a;
    ob = b;
    oc = c;
  }
  else if (std::string("bunge").compare(type) == 0) {
    oa = fmod(a - M_PI / 2.0, 2 * M_PI);
    ob = fmod(b, M_PI);
    oc = fmod(M_PI / 2.0 - c, 2 * M_PI);
  }
  else if (std::string("roe").compare(type) == 0) {
    oa = a;
    ob = b;
    oc = M_PI - c;
  }
  else {
    throw std::invalid_argument("Invalid Euler angle type, requires kocks, bunge, or roe.");
  }
}

void Orientation::from_kocks_(double a, double b, double c, double & oa,
                            double & ob, double & oc, std::string type)
{
  if (std::string("kocks").compare(type) == 0) {
    oa = a;
    ob = b;
    oc = c;
  }
  else if (std::string("bunge").compare(type) == 0) {
    oa = fmod(a + M_PI / 2.0, 2 * M_PI);
    ob = fmod(b, M_PI);
    oc = fmod(M_PI / 2.0 - c, 2 * M_PI);
  }
  else if (std::string("roe").compare(type) == 0) {
    oa = a;
    ob = b;
    oc = M_PI - c;
  }
  else {
    throw std::invalid_argument("Invalid Euler angle type, requires kocks, bunge, or roe.");
  }
}

void Orientation::kocks_to_matrix_(double a, double b, double c,
                                   double * const M)
{
  M[0] = -sin(c) * sin(a) - cos(c) * cos(a) * cos(b);
  M[1] = sin(c) * cos(a) - cos(c) * sin(a) * cos(b);
  M[2] = cos(c) * sin(b);
  M[3] = cos(c) * sin(a) - sin(c) * cos(a) * cos(b);
  M[4] = -cos(c) * cos(a) - sin(c) * sin(a) * cos(b);
  M[5] = sin(c) * sin(b);
  M[6] = cos(a) * sin(b);
  M[7] = sin(a) * sin(b);
  M[8] = cos(b);
}

Orientation operator*(const Orientation & lhs, const Orientation & rhs)
{
  Orientation cpy(lhs);
  cpy *= rhs;
  return cpy;
}

Orientation operator/(const Orientation & lhs, const Orientation & rhs)
{
  Orientation cpy(lhs);
  cpy /= rhs;
  return cpy;
}

Vector Orientation::apply(const Vector & a) const
{
  std::vector<double> qv(4);
  qv[0] = 0.0;
  qv[1] = a.data()[0];
  qv[2] = a.data()[1];
  qv[3] = a.data()[2];

  Quaternion q(qv);

  Quaternion base = Quaternion(std::vector<double>(quat_,quat_+4));
  Quaternion rv = base * q * base.conj();

  Vector res;
  std::copy(&rv.quat()[1], &rv.quat()[4], res.s());

  return res;
}

RankTwo Orientation::apply(const RankTwo & a) const
{
  RankTwo Q = to_tensor();
  return Q * (a * Q.transpose());
}

Symmetric Orientation::apply(const Symmetric & a) const
{
  Symmetric res;
  double * rd = res.s();
  const double * q = quat_;
  const double * d = a.data();

  double sq_two = sqrt(2.0);
  rd[0] = (-1 + 2*q[2]*q[2] + 2*q[3]*q[3])*(-1 + 2*q[2]*q[2] + 2*q[3]*q[3])*d[0] + 2*(2*q[1]*q[1]*(q[2]*q[2]*d[1] + q[3]*q[3]*d[2] + sq_two*q[2]*q[3]*d[3]) + q[1]*(2*q[0]*(2*q[2]*q[3]*(-d[1] + d[2]) + sq_two*q[2]*q[2]*d[3] - sq_two*q[3]*q[3]*d[3]) - sq_two*(-1 + 2*q[2]*q[2] + 2*q[3]*q[3])*(q[3]*d[4] + q[2]*d[5])) + q[0]*(2*q[0]*(q[3]*q[3]*d[1] + q[2]*q[2]*d[2] - sq_two*q[2]*q[3]*d[3]) - sq_two*(-1 + 2*q[2]*q[2] + 2*q[3]*q[3])*(q[2]*d[4] - q[3]*d[5])));
  rd[1] = (q[0]*q[1] - q[2]*q[3])*(4*(q[0]*q[1] - q[2]*q[3])*d[2] - sq_two*(d[3] - 2*q[1]*q[1]*d[3] - 2*q[3]*q[3]*d[3] + 2*q[1]*q[2]*d[4] + 2*q[0]*q[3]*d[4])) + (1 - 2*q[1]*q[1] - 2*q[3]*q[3])* ((1 - 2*q[1]*q[1] - 2*q[3]*q[3])*d[1] + sq_two*(-(q[0]*q[1]*d[3]) + q[2]*q[3]*d[3] + q[1]*q[2]*d[5] + q[0]*q[3]*d[5])) + (q[1]*q[2] + q[0]*q[3])*(4*(q[1]*q[2] + q[0]*q[3])*d[0] + sq_two*(d[5] - 2*(q[0]*q[1]*d[4] - q[2]*q[3]*d[4] + (q[1]*q[1] + q[3]*q[3])*d[5])));
  rd[2] = (1 - 2*q[1]*q[1] - 2*q[2]*q[2])* ((1 - 2*q[1]*q[1] - 2*q[2]*q[2])*d[2] + sq_two*(q[0]*q[1]*d[3] + q[2]*q[3]*d[3] - q[0]*q[2]*d[4] + q[1]*q[3]*d[4])) + (q[0]*q[1] + q[2]*q[3])*(4*(q[0]*q[1] + q[2]*q[3])*d[1] + sq_two*(d[3] - 2*q[1]*q[1]*d[3] - 2*q[2]*q[2]*d[3] - 2*q[0]*q[2]*d[5] + 2*q[1]*q[3]*d[5])) + (q[0]*q[2] - q[1]*q[3])*(4*(q[0]*q[2] - q[1]*q[3])*d[0] - sq_two*(d[4] - 2*q[1]*q[1]*d[4] - 2*q[2]*q[2]*d[4] + 2*q[0]*q[1]*d[5] + 2*q[2]*q[3]*d[5]));
  rd[3] = sq_two*(2*(-(q[0]*q[1]) + q[2]*q[3])* ((1 - 2*q[1]*q[1] - 2*q[2]*q[2])*d[2] + sq_two*(q[0]*q[1]*d[3] + q[2]*q[3]*d[3] - q[0]*q[2]*d[4] + q[1]*q[3]*d[4])) + ((-1 + 2*q[1]*q[1] + 2*q[3]*q[3])* (-4*(q[0]*q[1] + q[2]*q[3])*d[1] + sq_two*((-1 + 2*q[1]*q[1] + 2*q[2]*q[2])*d[3] + 2*(q[0]*q[2] - q[1]*q[3])*d[5]))) /2.0 + 2*(q[1]*q[2] + q[0]*q[3])* (-2*q[0]*q[2]*d[0] + 2*q[1]*q[3]*d[0] + d[4]/sq_two + sq_two*(-((q[1]*q[1] + q[2]*q[2])*d[4]) + (q[0]*q[1] + q[2]*q[3])*d[5])));
  rd[4] = sq_two*(2*(q[0]*q[2] + q[1]*q[3])* ((1 - 2*q[1]*q[1] - 2*q[2]*q[2])*d[2] + sq_two*(q[0]*q[1]*d[3] + q[2]*q[3]*d[3] - q[0]*q[2]*d[4] + q[1]*q[3]*d[4])) + 2*(q[1]*q[2] - q[0]*q[3])*(2*(q[0]*q[1] + q[2]*q[3])*d[1] + (d[3] - 2*q[1]*q[1]*d[3] - 2*q[2]*q[2]*d[3] - 2*q[0]*q[2]*d[5] + 2*q[1]*q[3]*d[5])/sq_two ) + ((-1 + 2*q[2]*q[2] + 2*q[3]*q[3])* (4*(q[0]*q[2] - q[1]*q[3])*d[0] + sq_two*((-1 + 2*q[1]*q[1] + 2*q[2]*q[2])*d[4] - 2*(q[0]*q[1] + q[2]*q[3])*d[5]))) /2.0);
  rd[5] = sq_two*(2*(q[0]*q[2] + q[1]*q[3])* (-2*q[0]*q[1]*d[2] + 2*q[2]*q[3]*d[2] + d[3]/sq_two + sq_two*(-((q[1]*q[1] + q[3]*q[3])*d[3]) + (q[1]*q[2] + q[0]*q[3])*d[4])) + 2*(q[1]*q[2] - q[0]*q[3])*((1 - 2*q[1]*q[1] - 2*q[3]*q[3])*d[1] + sq_two*(-(q[0]*q[1]*d[3]) + q[2]*q[3]*d[3] + q[1]*q[2]*d[5] + q[0]*q[3]*d[5])) + ((-1 + 2*q[2]*q[2] + 2*q[3]*q[3])* (-4*(q[1]*q[2] + q[0]*q[3])*d[0] + sq_two*(2*q[0]*q[1]*d[4] - 2*q[2]*q[3]*d[4] + (-1 + 2*q[1]*q[1] + 2*q[3]*q[3])*d[5])))/2.0);

  return Symmetric(apply(a.to_full()));
}

Skew Orientation::apply(const Skew & a) const
{
  Skew res;
  double * rw = res.s();
  const double * q = quat_;
  const double * w = a.data();

  rw[0] = -((1 - 2*std::pow(q[2],2) - 2*std::pow(q[3],2) + 4*std::pow(q[1],2)*(-1 + std::pow(q[0],2) + std::pow(q[1],2) + std::pow(q[2],2) + std::pow(q[3],2)))*-w[0] + 2*q[0]*(q[3]*w[1] + q[2]*-w[2]) - 2*q[1]*(-1 + 2*std::pow(q[0],2) + 2*std::pow(q[2],2) + 2*std::pow(q[3],2))* (q[2]*w[1] - q[3]*-w[2]) + std::pow(q[1],3)*(-4*q[2]*w[1] + 4*q[3]*-w[2]));
  rw[1] = (-4*(q[0]*q[2] + q[1]*q[3])*(q[0]*q[1]*-w[0] + q[2]*q[3]*-w[0] - q[0]*q[2]*w[1] + q[1]*q[3]*w[1]) - 2*(q[1]*q[2] - q[0]*q[3])*((-1 + 2*std::pow(q[1],2) + 2*std::pow(q[2],2))*-w[0] - 2*q[0]*q[2]*-w[2] + 2*q[1]*q[3]*-w[2]) + (1 - 2*std::pow(q[2],2) - 2*std::pow(q[3],2))* (w[1] - 2*std::pow(q[1],2)*w[1] - 2*std::pow(q[2],2)*w[1] + 2*q[0]*q[1]*-w[2] + 2*q[2]*q[3]*-w[2]));
  rw[2] = -(-2*q[0]*(q[2]*-w[0] + q[1]*w[1]) + 4*std::pow(q[0],2)*q[3]*(q[1]*-w[0] - q[2]*w[1] + q[3]*-w[2]) + (-1 + 2*std::pow(q[1],2) + 2*std::pow(q[2],2) + 2*std::pow(q[3],2))* (w[2] + 2*q[3]*(q[1]*-w[0] - q[2]*w[1] + q[3]*-w[2])));

  return res;
}

RankFour Orientation::apply(const RankFour & a) const
{
  RankTwo Q = to_tensor();

  RankFour res;
  for (size_t i = 0; i<3; i++) {
    for (size_t j = 0; j<3; j++) {
      for (size_t k = 0; k<3; k++) {
        for (size_t l = 0; l<3; l++) {
          for (size_t m = 0; m<3; m++) {
            for (size_t n = 0; n<3; n++) {
              for (size_t o = 0; o<3; o++) {
                for (size_t p = 0; p<3; p++) {
                  res(i,j,k,l) += Q(i,m) * Q(j,n) * Q(k,o) * Q(l,p) * a(m,n,o,p);
                }
              }
            }
          }
        }
      }
    }
  }
  return res;
}

SymSymR4 Orientation::apply(const SymSymR4 & a) const
{
  SymSymR4 res;

  double M[9];
  to_matrix(M);

  const double f1 = 1.0;
  const double f2 = sqrt(2.0);
  const double f3 = sqrt(2.0);
  double R[36] = {f1*M[0]*M[0],f1*M[1]*M[1],f1*M[2]*M[2],f3*M[1]*M[2],f3*M[2]*M[0],f3*M[0]*M[1],f1*M[3]*M[3],f1*M[4]*M[4],f1*M[5]*M[5],f3*M[4]*M[5],f3*M[5]*M[3],f3*M[3]*M[4],f1*M[6]*M[6],f1*M[7]*M[7],f1*M[8]*M[8],f3*M[7]*M[8],f3*M[8]*M[6],f3*M[6]*M[7],f2*M[3]*M[6],f2*M[4]*M[7],f2*M[5]*M[8],(M[4]*M[8]+M[5]*M[7]),(M[5]*M[6]+M[3]*M[8]),(M[3]*M[7]+M[4]*M[6]),f2*M[6]*M[0],f2*M[7]*M[1],f2*M[8]*M[2],(M[7]*M[2]+M[8]*M[1]),(M[8]*M[0]+M[6]*M[2]),(M[6]*M[1]+M[7]*M[0]),f2*M[0]*M[3],f2*M[1]*M[4],f2*M[2]*M[5],(M[1]*M[5]+M[2]*M[4]),(M[2]*M[3]+M[0]*M[5]),(M[0]*M[4]+M[1]*M[3])};

  rotate_matrix(6, 6, R, a.data(), res.s());

  return res;
}

double Orientation::distance(const Orientation & other) const
{
  double d = fabs(this->dot(other));
  if (d>1.0) d = 1.0;
  return acos(d);
}

std::vector<CrystalOrientation> random_orientations(int n)
{
  double u[3];
  double w, x, y, z;
  std::vector<CrystalOrientation> result;

  std::random_device rdev{};
  std::default_random_engine generator{rdev()};
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  for (int i=0; i<n; i++) {
    for (int j=0; j<3; j++) {
      u[j] = distribution(generator);
    }

    w = sqrt(1.0-u[0]) * sin(2.0 * M_PI * u[1]);
    x = sqrt(1.0-u[0]) * cos(2.0 * M_PI * u[1]);
    y = sqrt(u[0]) * sin(2.0 * M_PI * u[2]);
    z = sqrt(u[0]) * cos(2.0 * M_PI * u[2]);

    result.emplace_back(make_crystal_orientation(Orientation({w,x,y,z})));
  }

  return result;
}

Orientation wexp(const Skew & w)
{
  const double * wv = w.data();
  double x = norm2_vec(wv, 3);
  if (x == 0.0) {
    return Orientation();
  }
  double f = sin(x/2.0)/x;
  return Orientation({cos(x/2.0), f * wv[0], f * wv[1], f * wv[2]});
}

Skew wlog(const Orientation & q)
{
  Skew res;
  double s = q.quat()[0];
  const double * v = &q.quat()[1];

  std::copy(v, v+3, res.s());
  double f = 2.0 * acos(s / q.norm()) / norm2_vec(v, 3);

  return f * res;
}

double distance(const Orientation & q1, const Orientation & q2)
{
  return q1.distance(q2);
}

Orientation rotate_to(const Vector & a, const Vector & b)
{
  // Check that they are normalized
  Vector aa = a / a.norm();
  Vector bb = b / b.norm();

  // Get the axis
  Vector axis = aa.cross(bb).normalize();

  // Get the angle
  double angle = acos(aa.dot(bb));

  // Return the rotation
  return Orientation::createAxisAngle(axis.s(), angle);
}

Orientation rotate_to_family(const Vector & a, const Vector & b, double ang)
{
  Orientation base = rotate_to(a,b);

  Orientation null = Orientation::createAxisAngle(b.data(), ang);

  return null * base;
}

// Create from Euler angles stored in params
CrystalOrientation::CrystalOrientation(ParameterSet & params) :
    NEMLObject(params), Orientation()
{
  std::vector<double> angles = params.get_parameter<std::vector<double>>("angles");
  setEulerAngles(angles[0], angles[1], angles[2],
                 params.get_parameter<std::string>("angle_type"),
                 params.get_parameter<std::string>("angle_convention"));
}

std::string CrystalOrientation::type()
{
  return "CrystalOrientation";
}

ParameterSet CrystalOrientation::parameters()
{
  ParameterSet pset(CrystalOrientation::type());

  pset.add_parameter<std::vector<double>>("angles");
  pset.add_optional_parameter<std::string>("angle_type", std::string("radians"));
  pset.add_optional_parameter<std::string>("angle_convention",
                                           std::string("kocks"));

  return pset;
}

std::unique_ptr<NEMLObject> CrystalOrientation::initialize(ParameterSet & params)
{
  return neml::make_unique<CrystalOrientation>(params);
}

std::shared_ptr<CrystalOrientation> zero_orientation()
{
  ParameterSet params = CrystalOrientation::parameters();
  params.assign_parameter("angles", std::vector<double>({0,0,0}));

  return std::make_shared<CrystalOrientation>(params);
}

CrystalOrientation make_crystal_orientation(const Orientation & o)
{
  ParameterSet params = CrystalOrientation::parameters();
  double a, b, c;
  o.to_euler(a, b, c);
  params.assign_parameter("angles", std::vector<double>({a,b,c}));

  return CrystalOrientation(params);
}

} // namespace neml
