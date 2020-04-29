#ifndef DAMAGE_H
#define DAMAGE_H

#include "models.h"
#include "elasticity.h"
#include "larsonmiller.h"

#include "windows.h"

#include <memory>

namespace neml {

/// Small strain damage model
class NEML_EXPORT NEMLDamagedModel_sd: public NEMLModel_sd {
 public:
  /// Input is an elastic model, an undamaged base material, and the CTE
  NEMLDamagedModel_sd(
                      std::shared_ptr<LinearElasticModel> elastic,
                      std::shared_ptr<NEMLModel_sd> base,
                      std::shared_ptr<Interpolate> alpha,
                      bool truesdell);

  /// How many history variables?  Equal to base_history + ndamage
  virtual size_t nhist() const;
  /// Initialize base according to the base model and damage according to
  /// init_damage
  virtual int init_hist(double * const hist) const;

  /// The damaged stress update
  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n) = 0;

  /// Number of damage variables
  virtual size_t ndamage() const = 0;
  /// Setup the damage variables
  virtual int init_damage(double * const damage) const = 0;

  /// Override the elastic model
  virtual int set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);

 protected:
   std::shared_ptr<NEMLModel_sd> base_;
};

/// Scalar damage trial state
class SDTrialState: public TrialState {
 public:
  virtual ~SDTrialState() {};
  double e_np1[6];
  double e_n[6];
  double T_np1, T_n, t_np1, t_n, u_n, p_n;
  double s_n[6];
  double w_n;
  std::vector<double> h_n;
};

/// Special case where the damage variable is a scalar
class NEML_EXPORT NEMLScalarDamagedModel_sd: public NEMLDamagedModel_sd, public Solvable {
 public:
  /// Parameters are an elastic model, a base model, the CTE, a solver
  /// tolerance, the maximum number of solver iterations, and a verbosity
  /// flag
  NEMLScalarDamagedModel_sd(std::shared_ptr<LinearElasticModel> elastic,
                            std::shared_ptr<NEMLModel_sd> base,
                            std::shared_ptr<Interpolate> alpha,
                            double tol, int miter,
                            bool verbose, bool truesdell,
                            bool ekill, double dkill, double sfact);

  /// Stress update using the scalar damage model
  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  /// Equal to 1
  virtual size_t ndamage() const;
  /// Initialize to zero
  virtual int init_damage(double * const damage) const;

  /// Number of parameters for the solver
  virtual size_t nparams() const;
  /// Initialize the solver vector
  virtual int init_x(double * const x, TrialState * ts);
  /// The actual nonlinear residual and Jacobian to solve
  virtual int RJ(const double * const x, TrialState * ts,double * const R,
                 double * const J);
  /// Setup a trial state from known information
  int make_trial_state(const double * const e_np1, const double * const e_n,
                       double T_np1, double T_n, double t_np1, double t_n,
                       const double * const s_n, const double * const h_n,
                       double u_n, double p_n,
                       SDTrialState & tss);

  /// The scalar damage model
  virtual int damage(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const = 0;
  /// Derivative with respect to the damage variable
  virtual int ddamage_dd(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const = 0;
  /// Derivative with respect to the strain
  virtual int ddamage_de(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const = 0;
  /// Derivative with respect to the stress
  virtual int ddamage_ds(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const = 0;

 protected:
  int tangent_(const double * const e_np1, const double * const e_n,
               const double * const s_np1, const double * const s_n,
               double T_np1, double T_n, double t_np1, double t_n,
               double w_np1, double w_n, const double * const A_prime,
               double * const A);
  int ekill_update_(double T_np1, const double * const e_np1,
                    double * const s_np1,
                    double * const h_np1, const double * const h_n,
                    double * A_np1,
                    double & u_np1, double u_n,
                    double & p_np1, double p_n);

 protected:
  double tol_;
  int miter_;
  bool verbose_;
  bool ekill_;
  double dkill_;
  double sfact_;
};

/// Stack multiple scalar damage models together
class NEML_EXPORT CombinedDamageModel_sd: public NEMLScalarDamagedModel_sd {
 public:
  /// Parameters: elastic model, vector of damage models, the base model
  /// CTE, solver tolerance, solver max iterations, and a verbosity flag
  CombinedDamageModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::vector<std::shared_ptr<NEMLScalarDamagedModel_sd>> models,
      std::shared_ptr<NEMLModel_sd> base,
      std::shared_ptr<Interpolate> alpha,
      double tol, int miter,
      bool verbose, bool truesdell);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// The combined damage variable
  virtual int damage(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative with respect to damage
  virtual int ddamage_dd(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative with respect to strain
  virtual int ddamage_de(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative with respect to stress
  virtual int ddamage_ds(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;

  virtual int set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);


 protected:
  const std::vector<std::shared_ptr<NEMLScalarDamagedModel_sd>> models_;
};

static Register<CombinedDamageModel_sd> regCombinedDamageModel_sd;

/// Classical Hayhurst-Leckie-Rabotnov-Kachanov damage
class NEML_EXPORT ClassicalCreepDamageModel_sd: public NEMLScalarDamagedModel_sd {
 public:
  /// Parameters are the elastic model, the parameters A, xi, phi, the
  /// base model, the CTE, the solver tolerance, maximum iterations,
  /// and the verbosity flag.
  ClassicalCreepDamageModel_sd(
                            std::shared_ptr<LinearElasticModel> elastic,
                            std::shared_ptr<Interpolate> A,
                            std::shared_ptr<Interpolate> xi,
                            std::shared_ptr<Interpolate> phi,
                            std::shared_ptr<NEMLModel_sd> base,
                            std::shared_ptr<Interpolate> alpha,
                            double tol, int miter,
                            bool verbose, bool truesdell);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// The damage function d_np1 = d_n + (se / A)**xi (1 - d_np1)**(-phi) * dt
  virtual int damage(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt damage
  virtual int ddamage_dd(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt strain
  virtual int ddamage_de(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt stress
  virtual int ddamage_ds(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;

 protected:
  double se(const double * const s) const;

 protected:
  std::shared_ptr<Interpolate> A_;
  std::shared_ptr<Interpolate> xi_;
  std::shared_ptr<Interpolate> phi_;
};

static Register<ClassicalCreepDamageModel_sd> regClassicalCreepDamageModel_sd;

/// Base class of modular effective stresses used by ModularCreepDamageModel_sd
class NEML_EXPORT EffectiveStress: public NEMLObject {
 public:
  virtual int effective(const double * const s, double & eff) const = 0;
  virtual int deffective(const double * const s, double * const deff) const = 0;
};

/// von Mises stress
class NEML_EXPORT VonMisesEffectiveStress: public EffectiveStress
{
 public:
  VonMisesEffectiveStress();

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int effective(const double * const s, double & eff) const;
  virtual int deffective(const double * const s, double * const deff) const;
};

static Register<VonMisesEffectiveStress> regVonMisesEffectiveStress;

/// Mean stress
class NEML_EXPORT MeanEffectiveStress: public EffectiveStress
{
 public:
  MeanEffectiveStress();

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int effective(const double * const s, double & eff) const;
  virtual int deffective(const double * const s, double * const deff) const;
};

static Register<MeanEffectiveStress> regMeanEffectiveStress;

/// Huddleston stress
class NEML_EXPORT HuddlestonEffectiveStress: public EffectiveStress
{
 public:
  HuddlestonEffectiveStress(double b);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int effective(const double * const s, double & eff) const;
  virtual int deffective(const double * const s, double * const deff) const;

 private:
  double b_;
};

static Register<HuddlestonEffectiveStress> regHuddlestonEffectiveStress;

/// Maximum principal stress
class NEML_EXPORT MaxPrincipalEffectiveStress: public EffectiveStress
{
 public:
  MaxPrincipalEffectiveStress();

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int effective(const double * const s, double & eff) const;
  virtual int deffective(const double * const s, double * const deff) const;
};

static Register<MaxPrincipalEffectiveStress> regMaxPrincipalEffectiveStress;

/// Maximum of several effective stress measures
class NEML_EXPORT MaxSeveralEffectiveStress: public EffectiveStress
{
 public:
  MaxSeveralEffectiveStress(std::vector<std::shared_ptr<EffectiveStress>> measures);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int effective(const double * const s, double & eff) const;
  virtual int deffective(const double * const s, double * const deff) const;

 private:
  void select_(const double * const s, size_t & ind, double & value) const;

 private:
  std::vector<std::shared_ptr<EffectiveStress>> measures_;
};

static Register<MaxSeveralEffectiveStress> regMaxSeveralEffectiveStress;

/// Weighted some of several effective stresses
class NEML_EXPORT SumSeveralEffectiveStress: public EffectiveStress
{
 public:
  SumSeveralEffectiveStress(std::vector<std::shared_ptr<EffectiveStress>> measures,
                            std::vector<double> weights);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int effective(const double * const s, double & eff) const;
  virtual int deffective(const double * const s, double * const deff) const;

 private:
  std::vector<std::shared_ptr<EffectiveStress>> measures_;
  std::vector<double> weights_;
};

static Register<SumSeveralEffectiveStress> regSumSeveralEffectiveStress;

/// Modular version of Hayhurst-Leckie-Rabotnov-Kachanov damage
//    This model differs from the above in two ways
//      1) You can change the effective stress measure
//      2) There is an extra (1-w)^xi term in the formulation to make the
//         results match the old analytic solutions
class NEML_EXPORT ModularCreepDamageModel_sd: public NEMLScalarDamagedModel_sd {
 public:
  ModularCreepDamageModel_sd(
                            std::shared_ptr<LinearElasticModel> elastic,
                            std::shared_ptr<Interpolate> A,
                            std::shared_ptr<Interpolate> xi,
                            std::shared_ptr<Interpolate> phi,
                            std::shared_ptr<EffectiveStress> estress,
                            std::shared_ptr<NEMLModel_sd> base,
                            std::shared_ptr<Interpolate> alpha,
                            double tol, int miter,
                            bool verbose, bool truesdell,
                            bool ekill, double dkill,
                            double sfact);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// The damage function d_np1 = d_n + (se / A)**xi * (1-d_np1)**xi * (1 - d_np1)**(-phi) * dt
  virtual int damage(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt damage
  virtual int ddamage_dd(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt strain
  virtual int ddamage_de(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt stress
  virtual int ddamage_ds(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;

 protected:
  std::shared_ptr<Interpolate> A_;
  std::shared_ptr<Interpolate> xi_;
  std::shared_ptr<Interpolate> phi_;
  std::shared_ptr<EffectiveStress> estress_;
};

static Register<ModularCreepDamageModel_sd> regModularCreepDamageModel_sd;

/// Time-fraction ASME damage using a generic Larson-Miller relation and effective stress
class NEML_EXPORT LarsonMillerCreepDamageModel_sd: public NEMLScalarDamagedModel_sd {
 public:
  LarsonMillerCreepDamageModel_sd(
                            std::shared_ptr<LinearElasticModel> elastic,
                            std::shared_ptr<LarsonMillerRelation> lmr,
                            std::shared_ptr<EffectiveStress> estress,
                            std::shared_ptr<NEMLModel_sd> base,
                            std::shared_ptr<Interpolate> alpha,
                            double tol, int miter,
                            bool verbose, bool truesdell,
                            bool ekill, double dkill,
                            double sfact);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// The damage function d_np1 = d_n + 1/tr(s*(1-w), T) * dt
  virtual int damage(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt damage
  virtual int ddamage_dd(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt strain
  virtual int ddamage_de(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt stress
  virtual int ddamage_ds(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;

 protected:
  std::shared_ptr<LarsonMillerRelation> lmr_;
  std::shared_ptr<EffectiveStress> estress_;
};

static Register<LarsonMillerCreepDamageModel_sd> regLarsonMillerCreepDamageModel_sd;

/// A standard damage model where the damage rate goes as the plastic strain
class NEML_EXPORT NEMLStandardScalarDamagedModel_sd: public NEMLScalarDamagedModel_sd {
 public:
  /// Parameters: elastic model, base model, CTE, solver tolerance,
  /// solver maximum number of iterations, verbosity flag
  NEMLStandardScalarDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::shared_ptr<NEMLModel_sd> base,
      std::shared_ptr<Interpolate> alpha,
      double tol, int miter,
      bool verbose, bool truesdell);

  /// Damage, now only proportional to the inelastic effective strain
  virtual int damage(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt damage
  virtual int ddamage_dd(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt strain
  virtual int ddamage_de(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt stress
  virtual int ddamage_ds(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;

  /// The part of the damage rate proportional to the inelastic strain rate
  virtual int f(const double * const s_np1, double d_np1,
                double T_np1, double & f) const = 0;
  /// Derivative with respect to stress
  virtual int df_ds(const double * const s_np1, double d_np1,
                   double T_np1, double * const df) const = 0;
  /// Derivative with respect to damage
  virtual int df_dd(const double * const s_np1, double d_np1,
                   double T_np1, double & df) const = 0;

 protected:
  double dep(const double * const s_np1, const double * const s_n,
             const double * const e_np1, const double * const e_n,
             double T_np1) const;

};

/// Simple power law damage
class NEML_EXPORT NEMLPowerLawDamagedModel_sd: public NEMLStandardScalarDamagedModel_sd {
 public:
  /// Parameters are an elastic model, the constants A and a, the base
  /// material model, the CTE, a solver tolerance, solver maximum number
  /// of iterations, and a verbosity flag
  NEMLPowerLawDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::shared_ptr<Interpolate> A, std::shared_ptr<Interpolate> a,
      std::shared_ptr<NEMLModel_sd> base,
      std::shared_ptr<Interpolate> alpha,
      double tol, int miter,
      bool verbose, bool truesdell);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Damage = A * s_eq**a (times the inelastic strain rate)
  virtual int f(const double * const s_np1, double d_np1,
                double T_np1, double & f) const;
  /// Derivative of f wrt stress
  virtual int df_ds(const double * const s_np1, double d_np1, double T_np1,
                double * const df) const;
  /// Derivative of f wrt damage
  virtual int df_dd(const double * const s_np1, double d_np1, double T_np1,
                double & df) const;

 protected:
  double se(const double * const s) const;

 protected:
  std::shared_ptr<Interpolate> A_;
  std::shared_ptr<Interpolate> a_;
};

static Register<NEMLPowerLawDamagedModel_sd> regNEMLPowerLawDamagedModel_sd;

/// Simple exponential damage model
class NEML_EXPORT NEMLExponentialWorkDamagedModel_sd: public NEMLStandardScalarDamagedModel_sd {
 public:
  /// Parameters are the elastic model, parameters W0, k0, and af, the
  /// base material model, the CTE, a solver tolerance, maximum number
  /// of iterations, and a verbosity flag.
  NEMLExponentialWorkDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::shared_ptr<Interpolate> W0, std::shared_ptr<Interpolate> k0,
      std::shared_ptr<Interpolate> af,
      std::shared_ptr<NEMLModel_sd> base,
      std::shared_ptr<Interpolate> alpha,
      double tol, int miter,
      bool verbose, bool truesdell);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// damage rate is (d + k0)**af / W0 * s_eq
  virtual int f(const double * const s_np1, double d_np1,
                double T_np1, double & f) const;
  /// Derivative of damage wrt stress
  virtual int df_ds(const double * const s_np1, double d_np1, double T_np1,
                double * const df) const;
  /// Derivative of damage wrt damage
  virtual int df_dd(const double * const s_np1, double d_np1, double T_np1,
                double & df) const;

 protected:
  double se(const double * const s) const;

 protected:
  std::shared_ptr<Interpolate> W0_;
  std::shared_ptr<Interpolate> k0_;
  std::shared_ptr<Interpolate> af_;
};

static Register<NEMLExponentialWorkDamagedModel_sd> regNEMLExponentialWorkDamagedModel_sd;

class NEMLWorkRateFunctionDamage_sd: public NEMLScalarDamagedModel_sd {
 public:
  /// Parameters are the elastic model, the workfunction, parameters,
  /// the base model, the CTE, the solver tolerance, maximum iterations,
  /// and the verbosity flag.
  NEMLWorkRateFunctionDamage_sd(
                            std::shared_ptr<LinearElasticModel> elastic,
                            std::shared_ptr<Interpolate> workrate,
                            std::shared_ptr<Interpolate> Q,
                            std::shared_ptr<Interpolate> m,
                            std::shared_ptr<Interpolate> G,
                            std::shared_ptr<Interpolate> H,
                            std::shared_ptr<Interpolate> xi,
                            std::shared_ptr<Interpolate> phi,
                            std::shared_ptr<NEMLModel_sd> base,
                            std::shared_ptr<Interpolate> alpha,
                            double tol, int miter,
                            bool verbose, bool truesdell);

  /// String type for the object system
  static std::string type();
  /// Return the default parameters
  static ParameterSet parameters();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// damage d_np1 = d_n + dt*[pow(Wdot/Q,m)/f(Wdot) + G*pow(se/H,xi)/pow(1-d_np1,phi)]
  virtual int damage(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt damage
  virtual int ddamage_dd(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt strain
  virtual int ddamage_de(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  /// Derivative of damage wrt stress
  virtual int ddamage_ds(double d_np1, double d_n,
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;

 protected:
   double dep(const double * const s_np1, const double * const s_n,
              const double * const e_np1, const double * const e_n,
              double T_np1) const;

   double se(const double * const s) const;

 protected:
   std::shared_ptr<Interpolate> workrate_;
   std::shared_ptr<Interpolate> Q_;
   std::shared_ptr<Interpolate> m_;
   std::shared_ptr<Interpolate> G_;
   std::shared_ptr<Interpolate> H_;
   std::shared_ptr<Interpolate> xi_;
   std::shared_ptr<Interpolate> phi_;
};

static Register<NEMLWorkRateFunctionDamage_sd> regNEMLWorkRateFunctionDamage_sd;

} //namespace neml

#endif // DAMAGE_H
