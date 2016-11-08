#ifndef NEML_H
#define NEML_H

#include "elasticity.h"

#include <cstddef>
#include <memory>

namespace neml {

/// NEML material model interface definitions
//  All material models inherit from this base class.  It defines interfaces
//  and provides the methods for reading in material parameters.
class NEMLModel {
  public:
   NEMLModel();
   virtual ~NEMLModel();

   // To accommodate the three interfaces we need to store some
   // "secret" history variables
   virtual size_t nstore() const = 0;
   virtual int init_store(double * const store) const = 0;

   // These three interfaces are how FE programs can enter the model.
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const = 0;
   virtual int update_ldI(
       const double * const l_inc,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const = 0;
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const = 0;

   virtual size_t nhist() const = 0;  // Actual number of material history variables
   virtual int init_hist(double * const hist) const = 0;
};

/// Models implemented through the deformation gradient interface
//  
class NEMLModel_ldF: public NEMLModel {
  public:
   NEMLModel_ldF();
   virtual ~NEMLModel_ldF();
  
   // Up to the user to implement
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const = 0;

   // Defined here
   virtual size_t nstore() const;
   virtual int init_store(double * const store) const;
   virtual int update_ldI(
       const double * const l_inc,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const;
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const;

   // More interface
   virtual size_t nhist() const = 0;
   virtual int init_hist(double * const hist) const = 0;
};

/// Models implemented through the incremental large deformation interface
//
class NEMLModel_ldI: public NEMLModel {
  public:
   NEMLModel_ldI();
   virtual ~NEMLModel_ldI();
  
   // Up to the user to implement
   virtual int update_ldI(
       const double * const l_inc,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const = 0;

   // Defined here
   virtual size_t nstore() const;
   virtual int init_store(double * const store) const;
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const;
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const;

   // More interface
   virtual size_t nhist() const = 0;
   virtual int init_hist(double * const hist) const = 0;
};

/// Models implemented through the small deformation interface
//
class NEMLModel_sd: public NEMLModel {
  public:
    NEMLModel_sd();
    virtual ~NEMLModel_sd();

    // Up to the user to implement
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const = 0;

   // Defined here
   virtual size_t nstore() const;
   virtual int init_store(double * const store) const;
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const;
   virtual int update_ldI(
       const double * const l_inc,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const;

   // More interface
   virtual size_t nhist() const = 0;
   virtual int init_hist(double * const hist) const = 0;

};

/// Small strain linear elasticity as a test case
class SmallStrainElasticity: public NEMLModel_sd {
 public:
  SmallStrainElasticity(std::shared_ptr<LinearElasticModel> elastic);

   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) const;
   virtual size_t nhist() const;
   virtual int init_hist(double * const hist) const;

 private:
   std::shared_ptr<LinearElasticModel> elastic_;

};



} // namespace neml
#endif // NEML_H
