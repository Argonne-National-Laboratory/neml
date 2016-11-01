#ifndef DEVIATORIC_H
#define DEVIATORIC_H

#include "shear.h"
#include "surface.h"

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
  virtual int update(
      const double* const e_inc,
      double T_np1, double T_inc,
      double t_np1, double t_inc,
      double * const h_np1, const double * const h_n,
      double * const s_np1, const double * const s_n,
      double * const A_np1) const = 0;

};

/// Rate Independent Associative Flow Model
//  
//  This implements rate independent plasticity with a yield surface
//  defined as a f(deviatoric stress, history) and derivatives. 
//
//  This form lets me implement a whole bunch of rate independent models
//  quite quickly.  The algorithm used here is generalized closest point
//  projection.
//
//  See Simo & Hughes (1988) chapter 3 p. 146.
class RIAFModel: public DeviatoricModel {
 public:
  RIAFModel();
  virtual ~RIAFModel();

  // Defined here
  virtual size_t nhist() const;
  virtual int update(
      const double* const e_inc,
      double T_np1, double T_inc,
      double t_np1, double t_inc,
      double * const h_np1, const double * const h_n,
      double * const s_np1, const double * const s_n,
      double * const A_np1) const;

 private:
  std::unique_ptr<ShearModulus> shear_modulus_;
  std::unique_ptr<YieldSurface> yield_surface_;

};

} // namespace neml

#endif // DEVIATORIC_H
