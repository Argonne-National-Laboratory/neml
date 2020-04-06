#pragma once

#include "solvers.h"
#include "objects.h"
#include "elasticity.h"
#include "ri_flow.h"
#include "visco_flow.h"
#include "general_flow.h"
#include "interpolate.h"
#include "creep.h"

#include "windows.h"

#include <cstddef>
#include <memory>
#include <vector>
#include <cmath>
#include <iostream>

namespace neml {

/// NEML material model interface definitions
//  All material models inherit from this base class.  It defines interfaces
//  and provides the methods for reading in material parameters.
class NEML_EXPORT NEMLModel: public NEMLObject {
  public:
   /// Total number of stored internal variables
   virtual size_t nstore() const = 0;
   /// Initialize the internal variables
   virtual int init_store(double * const store) const = 0;

   /// Small strain update interface
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n) = 0;

   /// Large strain incremental update
   virtual int update_ld_inc(
       const double * const d_np1, const double * const d_n,
       const double * const w_np1, const double * const w_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1, double * const B_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n) = 0;

   /// Number of internal variables that are true material history
   virtual size_t nhist() const = 0;
   /// Initialize the history variables
   virtual int init_hist(double * const hist) const = 0;

   /// Instantaneous thermal expansion coefficient as a function of temperature
   virtual double alpha(double T) const = 0;
   /// Elastic strain for a given stress, temperature, and history state
   virtual int elastic_strains(const double * const s_np1,
                               double T_np1, const double * const h_np1,
                               double * const e_np1) const = 0;
};

/// Large deformation incremental update model
class NEML_EXPORT NEMLModel_ldi: public NEMLModel {
  public:
    NEMLModel_ldi();
    virtual ~NEMLModel_ldi();

   /// The small strain stress update interface
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n);

   /// Large strain incremental update
   virtual int update_ld_inc(
       const double * const d_np1, const double * const d_n,
       const double * const w_np1, const double * const w_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1, double * const B_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n) = 0;

   /// Number of stored variables
   virtual size_t nstore() const;
   /// Initialize stored variables
   virtual int init_store(double * const store) const;

   /// Number of stored variables that are true material history
   virtual size_t nhist() const = 0;
   /// Initialize the stored history
   virtual int init_hist(double * const hist) const = 0;
};

/// Small deformation stress update
class NEML_EXPORT NEMLModel_sd: public NEMLModel {
  public:
    /// All small strain models use small strain elasticity and CTE
    NEMLModel_sd(std::shared_ptr<LinearElasticModel> emodel,
                 std::shared_ptr<Interpolate> alpha,
                 bool truesdell);
    virtual ~NEMLModel_sd();

   /// The small strain stress update interface
   virtual int update_sd(
       const double * const e_np1, const double * const e_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n) = 0;

   /// Large strain incremental update
   virtual int update_ld_inc(
       const double * const d_np1, const double * const d_n,
       const double * const w_np1, const double * const w_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1, double * const B_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n);

   /// Number of stored variables
   virtual size_t nstore() const;
   /// Initialize stored variables
   virtual int init_store(double * const store) const;

   /// Number of stored variables that are true material history
   virtual size_t nhist() const = 0;
   /// Initialize the stored history
   virtual int init_hist(double * const hist) const = 0;

   /// Provide the instantaneous CTE
   virtual double alpha(double T) const;
   /// Returns the elasticity model, for sub-objects that want to use it
   const std::shared_ptr<const LinearElasticModel> elastic() const;

   /// Return the elastic strains
   virtual int elastic_strains(const double * const s_np1,
                               double T_np1, const double * const h_np1,
                               double * const e_np1) const;

   /// Used to override the linear elastic model to match another object's
   virtual int set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);

  private:
   int calc_tangent_(const double * const D, const double * const W,
                     const double * const C, const double * const S,
                     double * const A, double * const B);

  protected:
   std::shared_ptr<LinearElasticModel> elastic_;

  private:
   std::shared_ptr<Interpolate> alpha_;
   bool truesdell_;

};

/// Adaptive integration, tangent using the usual trick
class SubstepModel_sd: public NEMLModel_sd, public Solvable {
 public:
  SubstepModel_sd(std::shared_ptr<LinearElasticModel> emodel,
                  std::shared_ptr<Interpolate> alpha,
                  bool truesdell, double tol, int miter, bool verbose,
                  int max_divide, bool force_divide);

  /// Complete substep update
  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  /// Single step update
  virtual int update_step(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A, double * const E,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  /// Setup the trial state
  virtual TrialState * setup(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_n,
      const double * const h_n) = 0;

  /// Ignore update and take an elastic step
  virtual bool elastic_step(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_n,
      const double * const h_n) = 0;

  /// Interpret the x vector
  virtual int update_internal(
      const double * const x,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n) = 0;

  /// Minus the partial derivative of the residual with respect to the strain
  virtual int strain_partial(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_np1, const double * const s_n,
      const double * const h_np1, const double * const h_n,
      double * const de) = 0;

  /// Do the work calculation
  virtual int work_and_energy(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double & u_np1, double u_n,
      double & p_np1, double p_n) = 0;

 protected:
  double tol_;
  int miter_;
  bool verbose_;
  int max_divide_;
  bool force_divide_;
};

/// Small strain linear elasticity
//  This is generally only used as a basic test
class NEML_EXPORT SmallStrainElasticity: public NEMLModel_sd {
 public:
  /// Parameters are the minimum: an elastic model and a thermal expansion
  SmallStrainElasticity(std::shared_ptr<LinearElasticModel> elastic,
                        std::shared_ptr<Interpolate> alpha,
                        bool truesdell);

  /// Type for the object system
  static std::string type();
  /// Setup parameters for the object system
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Small strain stress update
  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n);
  /// Number of history variables (=0)
  virtual size_t nhist() const;
  /// Initialize history (none to setup)
  virtual int init_hist(double * const hist) const;
};

static Register<SmallStrainElasticity> regSmallStrainElasticity;

/// Small strain perfect plasticity trial state
//  Store data the solver needs and can be passed into solution interface
class SSPPTrialState : public TrialState {
 public:
  virtual ~SSPPTrialState() {};
  double ys, T;    // Yield stress, temperature
  double e_np1[6]; // Next strain
  double ep_n[6];  // Previous value of plastic strain
  double s_tr[6];  // Elastic trial stress
  double C[36];    // Stiffness
};

/// Small strain rate independent plasticity trial state
class SSRIPTrialState : public TrialState {
 public:
  virtual ~SSRIPTrialState() {};
  double ep_tr[6];          // Trial plastic strain
  double s_tr[6];           // Trial stress
  double e_np1[6];          // Next strain
  double C[36];             // Elastic stiffness
  double T;                 // Temperature
  std::vector<double> h_tr; // Trial history
};

/// Small strain creep+plasticity trial state
class SSCPTrialState : public TrialState {
 public:
  virtual ~SSCPTrialState() {};
  double ep_strain[6];            // Current plastic strain
  double e_n[6], e_np1[6];        // Previous and next total strain
  double s_n[6];                  // Previous stress
  double T_n, T_np1, t_n, t_np1;  // Next and previous time and temperature
  std::vector<double> h_n;        // Previous history vector
};

/// General inelastic integrator trial state
class GITrialState : public TrialState {
 public:
  virtual ~GITrialState() {};
  double e_dot[6];                // Strain rate
  double s_n[6];                  // Previous stress
  double T, Tdot, dt;             // Temperature, temperature rate, time inc.
  std::vector<double> h_n;        // Previous history
  double s_guess[6];              // Reasonable guess at the next stress
};

/// Small strain, associative, perfect plasticity
//    Algorithm is generalized closest point projection.
//    This degenerates to radial return for models where the gradient of
//    the yield surface is constant along lines from the origin to a point
//    in stress space outside the surface (i.e. J2).

class NEML_EXPORT SmallStrainPerfectPlasticity: public SubstepModel_sd {
 public:
  /// Parameters: elastic model, yield surface, yield stress, CTE,
  /// integration tolerance, maximum number of iterations,
  /// verbosity flag, and the maximum number of adaptive subdivisions
  SmallStrainPerfectPlasticity(std::shared_ptr<LinearElasticModel> elastic,
                               std::shared_ptr<YieldSurface> surface,
                               std::shared_ptr<Interpolate> ys,
                               std::shared_ptr<Interpolate> alpha,
                               double tol, int miter,
                               bool verbose,
                               int max_divide,
                               bool force_divide,
                               bool truesdell);

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Number of history variables (=0)
  virtual size_t nhist() const;
  /// Initialize history (nothing to do)
  virtual int init_hist(double * const hist) const;
  
  /// Number of nonlinear equations to solve in the integration
  virtual size_t nparams() const;
  /// Setup an initial guess for the nonlinear solution
  virtual int init_x(double * const x, TrialState * ts);
  /// Integration residual and jacobian equations
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

  /// Setup the trial state
  virtual TrialState * setup(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_n,
      const double * const h_n);
  
  /// Take an elastic step
  virtual bool elastic_step(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_n,
      const double * const h_n);

  /// Interpret the x vector
  virtual int update_internal(
      const double * const x,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n);

  /// Minus the partial derivative of the residual with respect to the strain
  virtual int strain_partial(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_np1, const double * const s_n,
      const double * const h_np1, const double * const h_n,
      double * de);

  /// Do the work calculation
  virtual int work_and_energy(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  /// Helper to return the yield stress
  double ys(double T) const;

  /// Setup a trial state for the solver from the input information
  int make_trial_state(const double * const e_np1, const double * const e_n,
                       double T_np1, double T_n, double t_np1, double t_n,
                       const double * const s_n, const double * const h_n,
                       SSPPTrialState & ts);

 private:
  std::shared_ptr<YieldSurface> surface_;
  std::shared_ptr<Interpolate> ys_;
};

static Register<SmallStrainPerfectPlasticity> regSmallStrainPerfectPlasticity;

/// Small strain, rate-independent plasticity
//    The algorithm used here is generalized closest point projection
//    for associative flow models.  For non-associative models the algorithm
//    may theoretically fail the discrete Kuhn-Tucker conditions, even
//    putting aside convergence issues on the nonlinear solver.
class NEML_EXPORT SmallStrainRateIndependentPlasticity: public SubstepModel_sd {
 public:
  /// Parameters: elasticity model, flow rule, CTE, solver tolerance, maximum
  /// solver iterations, verbosity flag, tolerance on the Kuhn-Tucker conditions
  /// check, and a flag on whether the KT conditions should be evaluated
  SmallStrainRateIndependentPlasticity(std::shared_ptr<LinearElasticModel> elastic,
                                       std::shared_ptr<RateIndependentFlowRule> flow,
                                       std::shared_ptr<Interpolate> alpha, bool truesdell,
                                       double tol, int miter, bool verbose,
                                       int max_divide, bool force_divide);

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Number of history variables
  virtual size_t nhist() const;
  /// Initialize history at time zero
  virtual int init_hist(double * const hist) const;

  /// Setup the trial state
  virtual TrialState * setup(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_n,
      const double * const h_n);

  /// Ignore update and take an elastic step
  virtual bool elastic_step(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_n,
      const double * const h_n);

  /// Interpret the x vector
  virtual int update_internal(
      const double * const x,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n);

  /// Minus the partial derivative of the residual with respect to the strain
  virtual int strain_partial(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_np1, const double * const s_n,
      const double * const h_np1, const double * const h_n,
      double * const de);

  /// Do the work calculation
  virtual int work_and_energy(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  /// Number of solver parameters
  virtual size_t nparams() const;
  /// Setup an iteration vector in the solver
  virtual int init_x(double * const x, TrialState * ts);
  /// Solver function returning the residual and jacobian of the nonlinear
  /// system of equations integrating the model
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

  /// Return the elastic model for subobjects
  const std::shared_ptr<const LinearElasticModel> elastic() const;

  /// Setup a trial state
  int make_trial_state(const double * const e_np1, const double * const e_n,
                       double T_np1, double T_n, double t_np1, double t_n,
                       const double * const s_n, const double * const h_n,
                       SSRIPTrialState & ts);

 private:
  std::shared_ptr<RateIndependentFlowRule> flow_;
};

static Register<SmallStrainRateIndependentPlasticity> regSmallStrainRateIndependentPlasticity;

/// Small strain, rate-independent plasticity + creep
//  Uses a combined iteration of a rate independent plastic + creep model
//  to solver overall update
class NEML_EXPORT SmallStrainCreepPlasticity: public NEMLModel_sd, public Solvable {
 public:
  /// Parameters are an elastic model, a base NEMLModel_sd, a CreepModel,
  /// the CTE, a solution tolerance, the maximum number of nonlinear
  /// iterations, a verbosity flag, and a scale factor to regularize
  /// the nonlinear equations.
  SmallStrainCreepPlasticity(
                             std::shared_ptr<LinearElasticModel> elastic,
                             std::shared_ptr<NEMLModel_sd> plastic,
                             std::shared_ptr<CreepModel> creep,
                             std::shared_ptr<Interpolate> alpha,
                             double tol, int miter,
                             bool verbose, double sf,
                             bool truesdell);

  /// Type for the object system
  static std::string type();
  /// Setup parameters for the object system
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Small strain stress update
  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  /// Number of history variables matches the base model
  virtual size_t nhist() const;
  /// Passes call for initial history to base model
  virtual int init_hist(double * const hist) const;

  /// The number of parameters in the nonlinear equation
  virtual size_t nparams() const;
  /// Initialize the nonlinear solver
  virtual int init_x(double * const x, TrialState * ts);
  /// Residual equation to solve and corresponding jacobian
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

  /// Setup a trial state from known information
  int make_trial_state(const double * const e_np1, const double * const e_n,
                       double T_np1, double T_n, double t_np1, double t_n,
                       const double * const s_n, const double * const h_n,
                       SSCPTrialState & ts);

  /// Set a new elastic model
  virtual int set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);

 private:
  int form_tangent_(double * const A, double * const B,
                    double * const A_np1);

 private:
  std::shared_ptr<NEMLModel_sd> plastic_;
  std::shared_ptr<CreepModel> creep_;

  double tol_, sf_;
  int miter_;
  bool verbose_;
};

static Register<SmallStrainCreepPlasticity> regSmallStrainCreepPlasticity;

/// Small strain general integrator
//    General NR one some stress rate + history evolution rate
//
class NEML_EXPORT GeneralIntegrator: public SubstepModel_sd {
 public:
  /// Parameters are an elastic model, a general flow rule,
  /// the CTE, the integration tolerance, the maximum
  /// nonlinear iterations, a verbosity flag, and the
  /// maximum number of subdivisions for adaptive integration
  GeneralIntegrator(std::shared_ptr<LinearElasticModel> elastic,
                    std::shared_ptr<GeneralFlowRule> rule,
                    std::shared_ptr<Interpolate> alpha,
                    bool truesdell, double tol, int miter, bool verbose,
                    int max_divide, bool force_divide);

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Setup the trial state
  virtual TrialState * setup(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_n,
      const double * const h_n);
  
  /// Take an elastic step
  virtual bool elastic_step(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_n,
      const double * const h_n);

  /// Interpret the x vector
  virtual int update_internal(
      const double * const x,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n);

  /// Minus the partial derivative of the residual with respect to the strain
  virtual int strain_partial(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      const double * const s_np1, const double * const s_n,
      const double * const h_np1, const double * const h_n,
      double * de);

  /// Do the work calculation
  virtual int work_and_energy(
      const TrialState * ts,
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  /// Number of history variables
  virtual size_t nhist() const;
  /// Initialize the history at time zero
  virtual int init_hist(double * const hist) const;

  /// Number of nonlinear equations
  virtual size_t nparams() const;
  /// Initialize a guess for the nonlinear iterations
  virtual int init_x(double * const x, TrialState * ts);
  /// The residual and jacobian for the nonlinear solve
  virtual int RJ(const double * const x, TrialState * ts,
                 double * const R, double * const J);

  /// Initialize a trial state
  int make_trial_state(const double * const e_np1, const double * const e_n,
                       double T_np1, double T_n, double t_np1, double t_n,
                       const double * const s_n, const double * const h_n,
                       GITrialState & ts);

  /// Set a new elastic model
  virtual int set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);

 private:
  std::shared_ptr<GeneralFlowRule> rule_;
};

static Register<GeneralIntegrator> regGeneralIntegrator;

/// Combines multiple small strain integrators based on regimes of
/// rate-dependent behavior.
//
//  This model uses the idea from Kocks & Mecking of a normalized activation
//  energy to call different integrators depending on the combination of
//  temperature and strain rate.
//
//  A typical use case would be switching from rate-independent to rate
//  dependent behavior based on a critical activation energy cutoff point
//
//  A user provides a vector of models (length n) and a corresponding vector
//  of normalized activation energies (length n-1) dividing the response into
//  segments.  All the models must have compatible hardening -- the history
//  is just going to be blindly passed between the models.
//
class NEML_EXPORT KMRegimeModel: public NEMLModel_sd {
 public:
  /// Parameters are an elastic model, a vector of valid NEMLModel_sd objects,
  /// the transition activation energies, the Boltzmann constant in appropriate
  /// units, a Burgers vector for normalization, a reference strain rate,
  /// and the CTE.
  KMRegimeModel(std::shared_ptr<LinearElasticModel> emodel,
                std::vector<std::shared_ptr<NEMLModel_sd>> models,
                std::vector<double> gs,
                double kboltz, double b, double eps0,
                std::shared_ptr<Interpolate> alpha,
                bool truesdell);

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// The small strain stress update
  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  /// The number of model history variables
  virtual size_t nhist() const;
  /// Initialize history at time zero
  virtual int init_hist(double * const hist) const;

  /// Set a new elastic model
  virtual int set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);

 private:
  double activation_energy_(const double * const e_np1,
                            const double * const e_n,
                            double T_np1,
                            double t_np1, double t_n);

 private:
  std::vector<std::shared_ptr<NEMLModel_sd>> models_;
  std::vector<double> gs_;
  double kboltz_, b_, eps0_;
};

static Register<KMRegimeModel> regKMRegimeModel;

} // namespace neml
