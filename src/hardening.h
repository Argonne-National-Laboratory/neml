#ifndef HARDENING_H
#define HARDENING_H

#include <cstddef>

namespace neml {

/// Interface for an associative hardening rule
//    Essentially a map between strain space internal variables and
//    stress space and its gradient.
//
//    The requirement of being associative is that the gradient (the
//    generalized plastic moduli D) is symmetric.
//
//    This isn't enforced here, but you should make sure your model
//    obeys it.
class AssociativeHardening {
 public:
  AssociativeHardening();
  virtual ~AssociativeHardening();
  
  // Interface
  virtual size_t nhist() const = 0;
  virtual int init_hist(double * const alpha) const = 0;
  virtual int q(const double * const alpha, double T, double * const qv) const = 0;
  virtual int D(const double * const alpha, double T, double * const Dv) const = 0;

  // Not required (can invert D) but nice if you can do it analytically
  virtual int D_inv(const double * const alpha, double T, double * const Dv) const;

};

/// IsoJ2 linear hardening
class IsoJ2LinearAHardening: public AssociativeHardening {
 public:
  IsoJ2LinearAHardening(double K0, double Kp);
  virtual ~IsoJ2LinearAHardening();

  virtual size_t nhist() const;
  virtual int init_hist(double * const alpha) const;
  virtual int q(const double * const alpha, double T, double * const qv) const;
  virtual int D(const double * const alpha, double T, double * const Dv) const;
  virtual int D_inv(const double * const alpha, double T, double * const Dv) const;

  // Getters
  double K0() const;
  double Kp() const;

 private:
  const double K0_, Kp_;

};



} // namespace neml

#endif // HARDENING_H
