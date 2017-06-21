#ifndef NEML_INTERFACE_H
#define NEML_INTERFACE_H

#include "nemlerror.h"

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
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n) = 0;
   virtual int update_ldI(
       const double * const l_inc,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n) = 0;
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n) = 0;

   virtual size_t nhist() const = 0;  // Actual number of material history variables
   virtual int init_hist(double * const hist) const = 0;

   // Well models will need this...
   virtual double alpha(double T) const = 0;
   // Helper for FEA output
   virtual int elastic_strains(const double * const s_np1,
                               double T_np1,
                               double * const e_np1) const = 0;

};

std::unique_ptr<NEMLModel> parse_xml(std::string fname, std::string mname, 
                                     int & ier);

} // namespace neml

#endif // NEML_INTERFACE_H
