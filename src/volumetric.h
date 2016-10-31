#ifndef VOLUMETRIC_H
#define VOLUMETRIC_H

#include <cstddef>

namespace neml {

/// Superclass for small strain models defining a volumetric response
class VolumetricModel {
 public:
  VolumetricModel();
  virtual ~VolumetricModel();

  // Wrapper to take 1 parameter (volumetric strain increment) form to 
  // full 6 vector form.  Note the oddity: this adds in the stress and 
  // tangent into the existing because it's supposed to update a deviator
  virtual int update(
      const double* const e_inc,
      double T_np1, double T_inc,
      double t_np1, double t_inc,
      double * const h_np1, const double * const h_n,
      double * const s_np1, const double * const s_n,
      double * const A_np1) const;

  // This is what the user needs to implement: the one parameter update
  virtual size_t nhist() const = 0;
  virtual int update_mean(
      double e_inc,
      double T_np1, double T_inc,
      double t_np1, double t_inc,
      double * const h_np1, const double * const h_n,
      double & s_np1, double s_n,
      double & A_np1) const = 0;
};

/// Basic volumetric model: constant bulk modulus
//    Parameters:
//      K     constant bulk modulus
//
//    nhist = 0
//
class VModelConstantK: public VolumetricModel {
 public:
  VModelConstantK(double K);
  virtual ~VModelConstantK();

  virtual size_t nhist() const;
  virtual int update_mean(
      double e_inc,
      double T_np1, double T_inc,
      double t_np1, double t_inc,
      double * const h_np1, const double * const h_n,
      double & s_np1, double s_n,
      double & A_np1) const;

  // Getters
  double K() const;

 private:
  const double K_;

};

} // namespace neml

#endif //VOLUMETRIC_H
