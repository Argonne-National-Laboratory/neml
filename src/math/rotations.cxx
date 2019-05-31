#include "rotations.h"

#include "nemlmath.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <chrono>

namespace neml {

Quaternion::Quaternion()
{
  quat_[0] = 1.0;
  quat_[1] = 0.0;
  quat_[2] = 0.0;
  quat_[3] = 0.0;
}

Quaternion::Quaternion(const std::vector<double> v)
{
  std::copy(&v[0], &v[4], quat_);
}

Quaternion::Quaternion(const double * const v)
{
  std::copy(v, v+4, quat_);
}

Quaternion::Quaternion(const Quaternion & other)
{
  std::copy(other.quat(), other.quat()+4, quat_);
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
  // Move
  std::copy(rhs.quat(), rhs.quat()+4, quat_);
  return *this;
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

// Unit quaternion stuff

Orientation Orientation::createVector(const double * const v)
{
  return Orientation(v);
}

Orientation Orientation::createRodrigues(const double * const r)
{
  double nr = 0.0;
  for (int i=0; i<3; i++) {
    nr += r[i] * r[i];
  }

  double f = 1.0 / (sqrt(1.0 + nr));
  
  double q[4];
  q[0] = f;
  q[1] = r[0] * f;
  q[2] = r[1] * f;
  q[3] = r[2] * f;

  return Orientation(q);
}

Orientation Orientation::createMatrix(const double * const M)
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
  
  double q[4];
  q[0] = s;
  q[1] = x;
  q[2] = y;
  q[3] = z;

  return Orientation(q);
}

Orientation Orientation::createAxisAngle(const double * const n, double a, 
                                         std::string angles)
{
  a = convert_angle(a, angles);
  double q[4];
  q[0] = cos(a/2.0);
  double s = sin(a/2.0);
  q[1] = s * n[0];
  q[2] = s * n[1];
  q[3] = s * n[2];

  return Orientation(q);
}

Orientation Orientation::createEulerAngles(double a, double b, double c,
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

  return Orientation::createMatrix(M);
}

// alpha in [0,2pi], beta in [0,pi] (sphere), gamma in [0,2pi] (circle)
Orientation Orientation::createHopf(double alpha, double beta, double gamma,
                                    std::string angles)
{
  alpha = convert_angle(alpha, angles);
  beta = convert_angle(beta, angles);
  gamma = convert_angle(gamma, angles);

  double q[4];
  q[0] = cos(beta / 2.0) * cos(gamma / 2.0);
  q[1] = cos(beta / 2.0) * sin(gamma / 2.0);
  q[2] = sin(beta / 2.0) * cos(alpha + gamma / 2.0);
  q[3] = sin(beta / 2.0) * sin(alpha + gamma / 2.0);

  return Orientation(q);
}

Orientation Orientation::createHyperspherical(double a1, double a2, double a3,
                                              std::string angles)
{
  a1 = convert_angle(a1, angles);
  a2 = convert_angle(a2, angles);
  a3 = convert_angle(a3, angles);

  double q[4];
  q[0] = cos(a1);
  q[1] = sin(a1)*cos(a2);
  q[2] = sin(a1)*sin(a2)*cos(a3);
  q[3] = sin(a1)*sin(a2)*sin(a3);

  return Orientation(q);
}

Orientation::Orientation() :
    Quaternion()
{

}

Orientation::Orientation(const double * const v) :
    Quaternion(v)
{
  normalize_();
}

Orientation::Orientation(const std::vector<double> v) :
    Quaternion(v)
{
  normalize_();
}

Orientation::Orientation(const Orientation & other) :
    Quaternion(other.quat())
{
  normalize_();
}

Orientation::Orientation(const Quaternion & other) :
    Quaternion(other.quat())
{
  normalize_();
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
  a = 2.0 * acos(quat_[0]);
  double s = sin(a / 2.0);
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
  double qv[4];
  qv[0] = 0.0;
  qv[1] = a.data()[0];
  qv[2] = a.data()[1];
  qv[3] = a.data()[2];

  Quaternion q(qv);
  Quaternion base = static_cast<Quaternion>(*this);
  Quaternion rv = base * q * base.conj();

  Vector res;
  std::copy(&rv.quat()[1], &rv.quat()[4], res.s());

  return res;
}

std::vector<Orientation> random_orientations(int n)
{
  double u[3];
  double w, x, y, z;
  std::vector<Orientation> result;
  
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

    result.emplace_back(Orientation({w,x,y,z}));
  }
  
  return result;
}

} // namespace cpfmwk
