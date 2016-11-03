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
  RIAFModel(ShearModulus & modulus, YieldSurface & surface, 
            AssociativeHardening & hardening);
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
  ShearModulus & modulus_;
  YieldSurface & surface_;
  AssociativeHardening & hardening_;

};

} // namespace neml

#endif // DEVIATORIC_H
