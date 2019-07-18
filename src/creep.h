#ifndef CREEP_H
#define CREEP_H

#include "objects.h"
#include "elasticity.h"
#include "interpolate.h"
#include "solvers.h"

namespace neml {

/// Scalar creep functions in terms of effective stress and strain
class ScalarCreepRule: public NEMLObject {
  public:
   /// Scalar creep strain rate as a function of effective stress, 
   /// effective strain, time, and temperature
   virtual int g(double seq, double eeq, double t, double T, double & g)
       const = 0;
   /// Derivative of scalar creep rate wrt effective stress
   virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) 
       const = 0;
   /// Derivative of scalar creep rate wrt effective strain
   virtual int dg_de(double seq, double eeq, double t, double T, double & dg) 
       const = 0;
   /// Derivative of scalar creep rate wrt time, defaults to zero
   virtual int dg_dt(double seq, double eeq, double t, double T, double & dg) 
       const;
   /// Derivative of scalar creep rate wrt temperature, defaults to zero
   virtual int dg_dT(double seq, double eeq, double t, double T, double & dg) 
       const;
};

/// Simple power law creep
class PowerLawCreep: public ScalarCreepRule {
 public:
  /// Parameters: prefector A and exponent n
  PowerLawCreep(std::shared_ptr<Interpolate> A, std::shared_ptr<Interpolate> n);
  
  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();
  
  /// rate = A * seq**n
  virtual int g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of rate wrt effective stress
  virtual int dg_ds(double seq, double eeq, double t, double T, double & dg)
      const;
  /// Derivative of rate wrt effective strain = 0
  virtual int dg_de(double seq, double eeq, double t, double T, double & dg)
      const;
  
  /// Getter for the prefactor
  double A(double T) const;
  /// Getter for the exponent
  double n(double T) const;

 private:
  const std::shared_ptr<const Interpolate> A_, n_;
};

static Register<PowerLawCreep> regPowerLawCreep;

/// A power law type model that uses KM concepts to switch between mechanisms
class RegionKMCreep: public ScalarCreepRule {
 public:
  /// Inputs: cuts in normalized activation energy, prefactors for each region
  /// exponents for each region, boltzmann constant, burgers vector, 
  /// reference strain rate, elastic model, to compute mu
  RegionKMCreep(std::vector<double> cuts, std::vector<double> A, 
                std::vector<double> B, double kboltz, double b, double eps0,
                std::shared_ptr<LinearElasticModel> emodel);
  
  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return the default parameters
  static ParameterSet parameters();
  
  /// See documentation for details of the creep rate
  virtual int g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of creep rate wrt effective stress
  virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const;
  /// Derivative of creep rate wrt effective strain
  virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const;

 private:
  void select_region_(double seq, double T, double & Ai, double & Bi) const;

 private:
  const std::vector<double> cuts_;
  const std::vector<double> A_;
  const std::vector<double> B_;
  const double kboltz_, b_, eps0_, b3_;
  const std::shared_ptr<LinearElasticModel> emodel_;
};

static Register<RegionKMCreep> regRegionKMCreep;

/// Classical Norton-Bailey creep
class NortonBaileyCreep: public ScalarCreepRule {
 public:
  /// Parameters: prefactor A, stress exponent n, time exponent m
  NortonBaileyCreep(std::shared_ptr<Interpolate> A, 
                    std::shared_ptr<Interpolate> m,
                    std::shared_ptr<Interpolate> n);
  
  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();
  
  /// creep rate = m * A**(1/m) * seq**(n/m) * eeq ** ((m-1)/m)
  virtual int g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of creep rate wrt effective stress
  virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const;
  /// Derivative of creep rate wrt effective strain
  virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const;
  
  /// Getter for the prefactor
  double A(double T) const;
  /// Getter for the stress exponent
  double m(double T) const;
  /// Getter for the time exponent
  double n(double T) const;

 private:
  const std::shared_ptr<const Interpolate> A_, m_, n_;
};

static Register<NortonBaileyCreep> regNortonBaileyCreep;

/// Classical Mukherjee creep
class MukherjeeCreep: public ScalarCreepRule {
 public:
  /// Parameters: elastic model (for shear modulus), prefactor, 
  /// stress exponent, reference lattice diffusivity, activation energy for
  /// lattice diffusion, burgers vector, boltzmann constant, gas constant
  MukherjeeCreep(std::shared_ptr<LinearElasticModel> emodel, double A, double n,
                 double D0, double Q, double b, double k, double R);
  
  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();
  
  /// scalar creep rate = A * D0 * exp(Q / (RT)) * mu * b / (k * T) * 
  /// (seq / mu)**n
  virtual int g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of creep rate wrt effective stress
  virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const;
  /// Derivative of creep rate wrt effective strain
  virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const;
  
  /// Getter for A
  double A() const;
  /// Getter for parameter n
  double n() const;
  /// Getter for parameter D0
  double D0() const;
  /// Getter for parameter Q
  double Q() const;
  /// Getter for parameter b
  double b() const;
  /// Getter for parameter k
  double k() const;
  /// Getter for parameter R
  double R() const;

 private:
  const std::shared_ptr<const LinearElasticModel> emodel_;
  const double A_, n_, D0_, Q_, b_, k_, R_;
};

static Register<MukherjeeCreep> regMukherjeeCreep;

/// A generic creep rate model where log(creep_rate) = Interpolate(log(stress))
class GenericCreep: public ScalarCreepRule {
 public:
  /// Parameters: interpolate giving the log creep rate
  GenericCreep(std::shared_ptr<Interpolate> cfn);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// scalar creep rate = exp(f(log(sigma))) 
  virtual int g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of creep rate wrt effective stress
  virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const;
  /// Derivative of creep rate wrt effective strain
  virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const;

 private:
  const std::shared_ptr<Interpolate> cfn_;
};

static Register<GenericCreep> regGenericCreep;

/// Creep trial state
class CreepModelTrialState : public TrialState {
 public:
  double T, dt, t;
  double s_np1[6];
  double e_n[6];
};

/// Interface to creep models
class CreepModel: public NEMLObject, public Solvable {
 public:
  /// Parameters are a solver tolerance, the maximum allowable iterations,
  /// and a verbosity flag
  CreepModel(double tol, int miter, bool verbose);
  
  /// Use the creep rate function to update the creep strain
  int update(const double * const s_np1, 
             double * const e_np1, const double * const e_n,
             double T_np1, double T_n,
             double t_np1, double t_n,
             double * const A_np1);
  
  /// The creep rate as a function of stress, strain, time, and temperature
  virtual int f(const double * const s, const double * const e, double t, double T, 
                double * const f) const = 0;
  /// The derivative of the creep rate wrt stress
  virtual int df_ds(const double * const s, const double * const e, double t, double T, 
                double * const df) const = 0;
  /// The derivative of the creep rate wrt strain
  virtual int df_de(const double * const s, const double * const e, double t, double T, 
                double * const df) const = 0;
  /// The derivative of the creep rate wrt time, defaults to zero
  virtual int df_dt(const double * const s, const double * const e, double t, double T, 
                double * const df) const;
  /// The derivative of the creep rate wrt temperature, defaults to zero
  virtual int df_dT(const double * const s, const double * const e, double t, double T, 
                double * const df) const;
  
  /// Setup a trial state for the solver
  int make_trial_state(const double * const s_np1, 
                       const double * const e_n,
                       double T_np1, double T_n,
                       double t_np1, double t_n,
                       CreepModelTrialState & ts) const;

  /// Number of solver parameters
  virtual size_t nparams() const;
  /// Setup the initial guess for the solver
  virtual int init_x(double * const x, TrialState * ts);
  /// The nonlinear residual and jacobian to solve
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

 private:
  int calc_tangent_(const double * const e_np1, CreepModelTrialState & ts, 
                    double * const A_np1);

 protected:
  const double tol_;
  const int miter_;
  const bool verbose_;
};

/// J2 creep based on a scalar creep rule
class J2CreepModel: public CreepModel {
 public:
  /// Parameters: scalar creep rule, nonlinear tolerance, maximum solver
  /// iterations, and a verbosity flag
  J2CreepModel(std::shared_ptr<ScalarCreepRule> rule,
               double tol, int miter, bool verbose);
  
  /// String type for the object system
  static std::string type();
  /// Initialize an object from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return the default parameters
  static ParameterSet parameters();
  
  /// creep_rate = dev(s) / ||dev(s)|| * scalar(effective_strain, 
  /// effective_stress, time, temperature)
  virtual int f(const double * const s, const double * const e, double t, double T, 
                double * const f) const;
  /// Derivative of creep rate wrt stress
  virtual int df_ds(const double * const s, const double * const e, double t, double T, 
                double * const df) const;
  /// Derivative of creep rate wrt strain
  virtual int df_de(const double * const s, const double * const e, double t, double T, 
                double * const df) const;
  /// Derivative of creep rate wrt time
  virtual int df_dt(const double * const s, const double * const e, double t, double T, 
                double * const df) const;
  /// Derivative of creep rate wrt temperature
  virtual int df_dT(const double * const s, const double * const e, double t, double T, 
                double * const df) const;

 private:
  // Helpers for computing the above
  double seq(const double * const s) const;
  double eeq(const double * const e) const;
  int sdir(double * const s) const;
  int edir(double * const e) const;

 private:
  std::shared_ptr<ScalarCreepRule> rule_;
};

static Register<J2CreepModel> regJ2CreepModel;

} // namespace neml

#endif // CREEP_H
