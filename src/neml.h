#ifndef NEML_H
#define NEML_H

#include <cstddef>
#include <memory>

#include "deviatoric.h"
#include "volumetric.h"

namespace neml {

/// NEML material model interface definitions
//  All material models inherit from this base class.  It defines interfaces
//  and provides the methods for reading in material parameters.
class NEMLModel {
  public:
   NEMLModel();
   virtual ~NEMLModel();

   virtual size_t nhist() = 0;  // Actual number of material history variables
   // To accommodate the three interfaces we need to store some 
   // "secret" history variables
   virtual size_t nstore() = 0;
  
   // These three interfaces are how FE programs can enter the model.
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1) = 0;
   virtual int update_ldI(
       const double* const l_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1) = 0;
   virtual int update_sd(
       const double* const e_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1)  = 0;
};

/// Models implemented through the deformation gradient interface
//  
class NEMLModel_ldF: public NEMLModel {
  public:
   NEMLModel_ldF();
   virtual ~NEMLModel_ldF();
  
   // Up to the user to implement
   virtual size_t nhist() = 0;
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1) = 0;

   // Defined here
   virtual size_t nstore();
   virtual int update_ldI(
       const double* const l_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1);
   virtual int update_sd(
       const double* const e_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1);
};

/// Models implemented through the incremental large deformation interface
//
class NEMLModel_ldI: public NEMLModel {
  public:
   NEMLModel_ldI();
   virtual ~NEMLModel_ldI();
  
   // Up to the user to implement
   virtual size_t nhist() = 0;
   virtual int update_ldI(
       const double* const l_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1) = 0;

   // Defined here
   virtual size_t nstore();
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1);
   virtual int update_sd(
       const double* const e_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1);
};

/// Models implemented through the small deformation interface
//
class NEMLModel_sd: public NEMLModel {
  public:
    NEMLModel_sd();
    virtual ~NEMLModel_sd();

    // Up to the user to implement
   virtual size_t nhist() = 0;
   virtual int update_sd(
       const double* const e_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1) = 0;

   // Defined here
   virtual size_t nstore();
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1);
   virtual int update_ldI(
       const double* const l_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1);

};

/// A small strain kinematics model with volumetric-deviatoric split
//  Uses two sub-objects: 
//    1) Deviatoric model
//    2) Volumetric model
class SplitModel_sd: public NEMLModel_sd {
  public:
   SplitModel_sd();
   virtual ~SplitModel_sd();

   virtual size_t nhist();
   virtual int update_sd(
       const double* const e_inc,
       double T_np1, double T_inc,
       double t_np1, double t_inc,
       double * const h_np1, const double * const h_n,
       double * const s_np1, const double * const s_n,
       double * const A_np1);

  private:
   std::unique_ptr<VolumetricModel> vol_model_;
   std::unique_ptr<DeviatoricModel> dev_model_;

};


} // namespace neml
#endif // NEML_H
