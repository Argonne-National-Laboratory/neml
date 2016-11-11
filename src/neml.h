#ifndef NEML_H
#define NEML_H

#include "elasticity.h"
#include "ri_flow.h"
#include "visco_flow.h"
#include "solvers.h"

#include <cstddef>
#include <memory>
#include <vector>

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
       double * const A_np1) = 0;
   virtual int update_ldI(
       const double * const l_inc,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) = 0;
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1) = 0;

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
       double * const A_np1) = 0;

   // Defined here
   virtual size_t nstore() const;
   virtual int init_store(double * const store) const;
   virtual int update_ldI(
       const double * const l_inc,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1);
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1);

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
       double * const A_np1) = 0;

   // Defined here
   virtual size_t nstore() const;
   virtual int init_store(double * const store) const;
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1);
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1);

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
       double * const A_np1) = 0;

   // Defined here
   virtual size_t nstore() const;
   virtual int init_store(double * const store) const;
   virtual int update_ldF(
       const double * const F_np1, const double * const F_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1);
   virtual int update_ldI(
       const double * const l_inc,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1);

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
      double * const A_np1);
  virtual size_t nhist() const;
  virtual int init_hist(double * const hist) const;

 private:
   std::shared_ptr<LinearElasticModel> elastic_;
};

/// Small strain, rate-independent plasticity
//    The algorithm used here is generalized closest point projection
//    for associative flow models.  For non-associative models the algorithm
//    may theoretically fail the discrete Kuhn-Tucker conditions, even
//    putting aside convergence issues on the nonlinear solver.
//
//    The class does check for Kuhn-Tucker violations when it returns, 
//    reporting an error if the conditions are violated.
class SmallStrainRateIndependentPlasticity: public NEMLModel_sd, public Solvable, public std::enable_shared_from_this<SmallStrainRateIndependentPlasticity> {
 public:
  SmallStrainRateIndependentPlasticity(std::shared_ptr<LinearElasticModel> elastic,
                                       std::shared_ptr<RateIndependentFlowRule> flow,
                                       double rtol = 1.0e-10, 
                                       double atol = 1.0e-12, int miter = 25,
                                       bool verbose = false, 
                                       double kttol = 1.0e-2,
                                       bool check_kt = false);
  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1);
  virtual size_t nhist() const;
  virtual int init_hist(double * const hist) const;

  virtual size_t nparams() const;
  virtual int init_x(double * const x);
  virtual int RJ(const double * const x, double * const R, double * const J);

  // Make this public for ease of testing
  int set_trial_state(const double * const e_np1, const double * const h_n,
                      double T_np1, double t_np1, double t_n);

 private:
  int calc_tangent_(const double * const x, const double * const s_np1,
                    const double * const h_np1, double dg, double * const A_np1);
  int check_K_T_(const double * const s_np1, const double * const h_np1, double T_np1, double dg);

  std::shared_ptr<LinearElasticModel> elastic_;
  std::shared_ptr<RateIndependentFlowRule> flow_;

  // Need to store the trial state in memory so that the solver 
  // can have access to it.  Not happy about it, but eh.
  double ep_tr_[6];
  double s_tr_[6];
  double e_np1_[6];
  double C_[36];
  double T_;
  std::vector<double> h_tr_;
  double rtol_, atol_, kttol_;
  int miter_;
  bool verbose_, check_kt_;
};

/// Small strain viscoplasticity
//    The NR algorithm here is a closest point projection.
//
//    Note this is very nearly the same algorithm as above.
//

class SmallStrainViscoPlasticity: public NEMLModel_sd, public Solvable, public std::enable_shared_from_this<SmallStrainViscoPlasticity> {
 public:
  SmallStrainViscoPlasticity(std::shared_ptr<LinearElasticModel> elastic,
                             std::shared_ptr<ViscoPlasticFlowRule> flow,
                             double rtol = 1.0e-6, 
                             double atol = 1.0e-10, int miter = 250,
                             bool verbose = false);
  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1);
  virtual size_t nhist() const;
  virtual int init_hist(double * const hist) const;

  virtual size_t nparams() const;
  virtual int init_x(double * const x);
  virtual int RJ(const double * const x, double * const R, double * const J);

  // Make this public for ease of testing
  int set_trial_state(const double * const e_np1, const double * const h_n,
                      double T_np1, double t_np1, double t_n);

 private:
  int calc_tangent_(const double * const x, const double * const s_np1,
                    const double * const h_np1, double dg, double * const A_np1);

  std::shared_ptr<LinearElasticModel> elastic_;
  std::shared_ptr<ViscoPlasticFlowRule> flow_;

  // Need to store the trial state in memory so that the solver 
  // can have access to it.  Not happy about it, but eh.
  double ep_tr_[6];
  double s_tr_[6];
  double e_np1_[6];
  double C_[36];
  double T_, dt_;
  std::vector<double> h_tr_;
  double rtol_, atol_;
  int miter_;
  bool verbose_;
};



} // namespace neml
#endif // NEML_H
