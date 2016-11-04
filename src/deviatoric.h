#ifndef DEVIATORIC_H
#define DEVIATORIC_H

#include "shear.h"
#include "surface.h"
#include "hardening.h"

#include <cstddef>
#include <memory>

namespace neml {

/// Superclass for small strain models defining a deviatoric response
class DeviatoricModel {
 public:
  DeviatoricModel();
  virtual ~DeviatoricModel();

  // Up to the user to implement
  virtual size_t nhist() const = 0;
  virtual int init_hist(double* const h) const = 0;
  virtual int update(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1) const = 0;

 protected:
  int get_dev(const double * const in, double * const out) const;

};

/// Linear elastic model
//
//  Very straightforward
class LEModel: public DeviatoricModel {
 public:
  LEModel(std::shared_ptr<ShearModulus> modulus);
  virtual ~LEModel();

  // Defined interface
  virtual size_t nhist() const;
  virtual int init_hist(double* const h) const;
  virtual int update(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1) const;

 private:
  std::shared_ptr<ShearModulus> modulus_;

};

/// Rate Independent Associative Flow Model
//  
//  This implements rate independent plasticity using surface and 
//  associative hardening objects.
//
//  This form lets me implement a whole bunch of rate independent models
//  quite quickly.  The algorithm used here is generalized closest point
//  projection.
//
//  See Simo & Hughes (1988) chapter 3 p. 146.
class RIAFModel: public DeviatoricModel {
 public:
  RIAFModel(std::shared_ptr<ShearModulus> modulus, 
            std::shared_ptr<YieldSurface> surface, 
            std::shared_ptr<AssociativeHardening> hardening,
            double tol = 1.0e-6,
            int miter = 25);
  virtual ~RIAFModel();

  // Defined here
  virtual size_t nhist() const;
  virtual int init_hist(double* const h) const;
  virtual int update(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1) const;

 private:
  bool take_step(const double * const e_dev, double mu, double T,
                 double * const ep_np1, double * const alpha_np1,
                 const double * const ep_n, const double * const alpha_n,
                 double & dg, double * const s) const;
  void get_jacobian(const double * const s, const double * const q,
                    const double * const Di, double T,
                    double mu, double dg, double * const J) const;

  std::shared_ptr<ShearModulus> modulus_;
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<AssociativeHardening> hardening_;

  const double tol_;
  const int miter_;
  bool verbose_;

};

} // namespace neml

#endif // DEVIATORIC_H
