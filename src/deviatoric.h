#ifndef DEVIATORIC_H
#define DEVIATORIC_H

#include <cstddef>

namespace neml {

/// Superclass for small strain models defining a deviatoric response
class DeviatoricModel {
 public:
  DeviatoricModel();
  virtual ~DeviatoricModel();

  // Up to the user to implement
  virtual size_t nhist() = 0;
  virtual int update(
      const double* const e_inc,
      double T_np1, double T_inc,
      double t_np1, double t_inc,
      double * const h_np1, const double * const h_n,
      double * const s_np1, const double * const s_n,
      double * const A_np1) = 0;

};

} // namespace neml

#endif // DEVIATORIC_H
