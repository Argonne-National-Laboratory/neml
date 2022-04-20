#ifndef CREEP_H
#define CREEP_H

#include "objects.h"
#include "elasticity.h"
#include "interpolate.h"
#include "solvers.h"

#include "windows.h"

namespace neml {

/// Scalar creep functions in terms of effective stress and strain
class NEML_EXPORT ScalarCreepRule: public NEMLObject {
  public:
   ScalarCreepRule(ParameterSet & params);
   /// Scalar creep strain rate as a function of effective stress,
   /// effective strain, time, and temperature
   virtual void g(double seq, double eeq, double t, double T, double & g)
       const = 0;
   /// Derivative of scalar creep rate wrt effective stress
   virtual void dg_ds(double seq, double eeq, double t, double T, double & dg)
       const = 0;
   /// Derivative of scalar creep rate wrt effective strain
   virtual void dg_de(double seq, double eeq, double t, double T, double & dg)
       const = 0;
   /// Derivative of scalar creep rate wrt time, defaults to zero
   virtual void dg_dt(double seq, double eeq, double t, double T, double & dg)
       const;
   /// Derivative of scalar creep rate wrt temperature, defaults to zero
   virtual void dg_dT(double seq, double eeq, double t, double T, double & dg)
       const;
};

/// Creep rate law from Blackburn 1972
class NEML_EXPORT BlackburnMinimumCreep: public ScalarCreepRule {
 public:
  BlackburnMinimumCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// rate = A * (sinh(beta*s/n)^n * exp(-Q/(R*T))
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg)
      const;
  /// Derivative of rate wrt effective strain = 0
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg)
      const;
  /// Derivative of rate wrt temperature
  virtual void dg_dT(double seq, double eeq, double t, double T, double & dg)
      const;

 private:
  const std::shared_ptr<const Interpolate> A_, n_, beta_;
  const double R_, Q_;
};

static Register<BlackburnMinimumCreep> regBlackburnMinimumCreep;

/// Creep rate law from Swindeman 1999
class NEML_EXPORT SwindemanMinimumCreep: public ScalarCreepRule {
 public:
  SwindemanMinimumCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// rate = C * S^n * exp(V*S) * exp(-Q/T)
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg)
      const;
  /// Derivative of rate wrt effective strain = 0
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg)
      const;
  /// Derivative of rate wrt temperature
  virtual void dg_dT(double seq, double eeq, double t, double T, double & dg)
      const;

 private:
  const double C_, n_, V_, Q_, shift_;
};

static Register<SwindemanMinimumCreep> regSwindemanMinimumCreep;

/// Simple power law creep
class NEML_EXPORT PowerLawCreep: public ScalarCreepRule {
 public:
  /// Parameters: prefector A and exponent n
  PowerLawCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// rate = A * seq**n
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg)
      const;
  /// Derivative of rate wrt effective strain = 0
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg)
      const;

  /// Getter for the prefactor
  double A(double T) const;
  /// Getter for the exponent
  double n(double T) const;

 private:
  const std::shared_ptr<const Interpolate> A_, n_;
};

static Register<PowerLawCreep> regPowerLawCreep;

/// Power law creep normalized in a slightly nicer way
class NEML_EXPORT NormalizedPowerLawCreep: public ScalarCreepRule {
 public:
  /// Parameters: stress divisor s0 and exponent n
  NormalizedPowerLawCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// rate = (seq/s0)**n
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg)
      const;
  /// Derivative of rate wrt effective strain = 0
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg)
      const;

 private:
  const std::shared_ptr<const Interpolate> s0_, n_;
};

static Register<NormalizedPowerLawCreep> regNormalizedPowerLawCreep;

/// A power law type model that uses KM concepts to switch between mechanisms
class NEML_EXPORT RegionKMCreep: public ScalarCreepRule {
 public:
  /// Inputs: cuts in normalized activation energy, prefactors for each region
  /// exponents for each region, boltzmann constant, burgers vector,
  /// reference strain rate, elastic model, to compute mu
  RegionKMCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return the default parameters
  static ParameterSet parameters();

  /// See documentation for details of the creep rate
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of creep rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg) const;
  /// Derivative of creep rate wrt effective strain
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg) const;

 private:
  void select_region_(double seq, double T, double & Ai, double & Bi) const;

 private:
  const std::vector<double> cuts_;
  const std::vector<std::shared_ptr<Interpolate>> A_;
  const std::vector<std::shared_ptr<Interpolate>> B_;
  const double kboltz_, b_, eps0_, b3_;
  const std::shared_ptr<LinearElasticModel> emodel_;
  double shift_;
};

static Register<RegionKMCreep> regRegionKMCreep;

/// Classical Norton-Bailey creep
class NEML_EXPORT NortonBaileyCreep: public ScalarCreepRule {
 public:
  /// Parameters: prefactor A, stress exponent n, time exponent m
  NortonBaileyCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// creep rate = m * A**(1/m) * seq**(n/m) * eeq ** ((m-1)/m)
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of creep rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg) const;
  /// Derivative of creep rate wrt effective strain
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg) const;

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
class NEML_EXPORT MukherjeeCreep: public ScalarCreepRule {
 public:
  /// Parameters: elastic model (for shear modulus), prefactor,
  /// stress exponent, reference lattice diffusivity, activation energy for
  /// lattice diffusion, burgers vector, boltzmann constant, gas constant
  MukherjeeCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// scalar creep rate = A * D0 * exp(Q / (RT)) * mu * b / (k * T) *
  /// (seq / mu)**n
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of creep rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg) const;
  /// Derivative of creep rate wrt effective strain
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg) const;

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
class NEML_EXPORT GenericCreep: public ScalarCreepRule {
 public:
  /// Parameters: interpolate giving the log creep rate
  GenericCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// scalar creep rate = exp(f(log(sigma)))
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of creep rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg) const;
  /// Derivative of creep rate wrt effective strain
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg) const;

 private:
  const std::shared_ptr<Interpolate> cfn_;
};

static Register<GenericCreep> regGenericCreep;

/// The hopelessly complicated 2.25Cr minimum creep rate model
//  Units: MPa and hours
class NEML_EXPORT MinCreep225Cr1MoCreep: public ScalarCreepRule {
 public:
  /// Parameters: prefector A and exponent n
  MinCreep225Cr1MoCreep(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Setup from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// See documentation for the hideous formula
  virtual void g(double seq, double eeq, double t, double T, double & g) const;
  /// Derivative of rate wrt effective stress
  virtual void dg_ds(double seq, double eeq, double t, double T, double & dg)
      const;
  /// Derivative of rate wrt effective strain = 0
  virtual void dg_de(double seq, double eeq, double t, double T, double & dg)
      const;

 private:
  double e1_(double seq, double T) const;
  double e2_(double seq, double T) const;

  double de1_(double seq, double T) const;
  double de2_(double seq, double T) const;

 private:
  const static std::shared_ptr<PiecewiseLinearInterpolate> U;
};

static Register<MinCreep225Cr1MoCreep> regMinCreep225Cr1MoCreep;

/// Creep trial state
class CreepModelTrialState : public TrialState {
 public:
  virtual ~CreepModelTrialState() {};
  double T, dt, t;
  double s_np1[6];
  double e_n[6];
};

/// Interface to creep models
class NEML_EXPORT CreepModel: public NEMLObject, public Solvable {
 public:
  /// Parameters are a solver tolerance, the maximum allowable iterations,
  /// and a verbosity flag
  CreepModel(ParameterSet & params);

  /// Use the creep rate function to update the creep strain
  void update(const double * const s_np1,
             double * const e_np1, const double * const e_n,
             double T_np1, double T_n,
             double t_np1, double t_n,
             double * const A_np1);

  /// The creep rate as a function of stress, strain, time, and temperature
  virtual void f(const double * const s, const double * const e, double t, double T,
                double * const f) const = 0;
  /// The derivative of the creep rate wrt stress
  virtual void df_ds(const double * const s, const double * const e, double t, double T,
                double * const df) const = 0;
  /// The derivative of the creep rate wrt strain
  virtual void df_de(const double * const s, const double * const e, double t, double T,
                double * const df) const = 0;
  /// The derivative of the creep rate wrt time, defaults to zero
  virtual void df_dt(const double * const s, const double * const e, double t, double T,
                double * const df) const;
  /// The derivative of the creep rate wrt temperature, defaults to zero
  virtual void df_dT(const double * const s, const double * const e, double t, double T,
                double * const df) const;

  /// Setup a trial state for the solver
  void make_trial_state(const double * const s_np1,
                       const double * const e_n,
                       double T_np1, double T_n,
                       double t_np1, double t_n,
                       CreepModelTrialState & ts) const;

  /// Number of solver parameters
  virtual size_t nparams() const;
  /// Setup the initial guess for the solver
  virtual void init_x(double * const x, TrialState * ts);
  /// The nonlinear residual and jacobian to solve
  virtual void RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

 private:
  void calc_tangent_(const double * const e_np1, CreepModelTrialState & ts,
                    double * const A_np1);

 protected:
  const double rtol_, atol_;
  const int miter_;
  const bool verbose_, linesearch_;
};

/// J2 creep based on a scalar creep rule
class NEML_EXPORT J2CreepModel: public CreepModel {
 public:
  /// Parameters: scalar creep rule, nonlinear tolerance, maximum solver
  /// iterations, and a verbosity flag
  J2CreepModel(ParameterSet & params);

  /// String type for the object system
  static std::string type();
  /// Initialize an object from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return the default parameters
  static ParameterSet parameters();

  /// creep_rate = dev(s) / ||dev(s)|| * scalar(effective_strain,
  /// effective_stress, time, temperature)
  virtual void f(const double * const s, const double * const e, double t, double T,
                double * const f) const;
  /// Derivative of creep rate wrt stress
  virtual void df_ds(const double * const s, const double * const e, double t, double T,
                double * const df) const;
  /// Derivative of creep rate wrt strain
  virtual void df_de(const double * const s, const double * const e, double t, double T,
                double * const df) const;
  /// Derivative of creep rate wrt time
  virtual void df_dt(const double * const s, const double * const e, double t, double T,
                double * const df) const;
  /// Derivative of creep rate wrt temperature
  virtual void df_dT(const double * const s, const double * const e, double t, double T,
                double * const df) const;

 private:
  // Helpers for computing the above
  double seq(const double * const s) const;
  double eeq(const double * const e) const;
  void sdir(double * const s) const;
  void edir(double * const e) const;

 private:
  std::shared_ptr<ScalarCreepRule> rule_;
};

static Register<J2CreepModel> regJ2CreepModel;

} // namespace neml

#endif // CREEP_H
