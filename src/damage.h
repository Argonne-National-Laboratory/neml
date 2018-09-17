#ifndef DAMAGE_H
#define DAMAGE_H

#include "models.h"
#include "elasticity.h"

#include <memory>

namespace neml {

/// Small strain damage model
class NEMLDamagedModel_sd: public NEMLModel_sd {
 public:
  NEMLDamagedModel_sd(
                      std::shared_ptr<LinearElasticModel> elastic,
                      std::shared_ptr<NEMLModel_sd> base,
                      std::shared_ptr<Interpolate> alpha = nullptr);
  NEMLDamagedModel_sd(std::shared_ptr<LinearElasticModel> elastic,
                      std::shared_ptr<NEMLModel_sd> base,
                      double alpha = 0.0);

  virtual size_t nhist() const;
  virtual int init_hist(double * const hist) const;

  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n) = 0;

  virtual size_t ndamage() const = 0;
  virtual int init_damage(double * const damage) const = 0;

  virtual int set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);

 protected:
   std::shared_ptr<NEMLModel_sd> base_;

};

/// Scalar damage trial state
class SDTrialState: public TrialState {
 public:
  double e_np1[6];
  double e_n[6];
  double T_np1, T_n, t_np1, t_n, u_n, p_n;
  double s_n[6];
  double w_n;
  std::vector<double> h_n;
};

/// Special case where the damage variable is a scalar
class NEMLScalarDamagedModel_sd: public NEMLDamagedModel_sd, public Solvable {
 public:
  NEMLScalarDamagedModel_sd(std::shared_ptr<LinearElasticModel> elastic,
                            std::shared_ptr<NEMLModel_sd> base,
                            std::shared_ptr<Interpolate> alpha = nullptr,
                            double tol = 1.0e-8, int miter = 50,
                            bool verbose = false);
  NEMLScalarDamagedModel_sd(std::shared_ptr<LinearElasticModel> elastic,
                            std::shared_ptr<NEMLModel_sd> base,
                            double alpha = 0.0,
                            double tol = 1.0e-8, int miter = 50,
                            bool verbose = false);

  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n);

  virtual size_t ndamage() const;
  virtual int init_damage(double * const damage) const;

  virtual size_t nparams() const;
  virtual int init_x(double * const x, TrialState * ts);
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);
  int make_trial_state(const double * const e_np1, const double * const e_n,
                       double T_np1, double T_n, double t_np1, double t_n,
                       const double * const s_n, const double * const h_n,
                       double u_n, double p_n,
                       SDTrialState & tss);

  virtual int damage(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const = 0;
  virtual int ddamage_dd(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const = 0;
  virtual int ddamage_de(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const = 0;
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

 protected:
  double tol_;
  int miter_;
  bool verbose_;
};

class CombinedDamageModel_sd: public NEMLScalarDamagedModel_sd {
 public:
  CombinedDamageModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::vector<std::shared_ptr<NEMLScalarDamagedModel_sd>> models,
      std::shared_ptr<NEMLModel_sd> base,
      std::shared_ptr<Interpolate> alpha = nullptr,
      double tol = 1.0e-8, int miter = 50,
      bool verbose = false);
  CombinedDamageModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::vector<std::shared_ptr<NEMLScalarDamagedModel_sd>> models,
      std::shared_ptr<NEMLModel_sd> base,
      double alpha = 0.0,
      double tol = 1.0e-8, int miter = 50,
      bool verbose = false);

  static std::string type();
  static ParameterSet parameters();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int damage(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  virtual int ddamage_dd(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  virtual int ddamage_de(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
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

class ClassicalCreepDamageModel_sd: public NEMLScalarDamagedModel_sd {
 public:
  ClassicalCreepDamageModel_sd(
                            std::shared_ptr<LinearElasticModel> elastic,
                            std::shared_ptr<Interpolate> A,
                            std::shared_ptr<Interpolate> xi,
                            std::shared_ptr<Interpolate> phi,
                            std::shared_ptr<NEMLModel_sd> base,
                            std::shared_ptr<Interpolate> alpha = nullptr,
                            double tol = 1.0e-8, int miter = 50,
                            bool verbose = false);
  ClassicalCreepDamageModel_sd(
                            std::shared_ptr<LinearElasticModel> elastic,
                            double A, double xi, double phi,
                            std::shared_ptr<NEMLModel_sd> base,
                            double alpha = 0.0,
                            double tol = 1.0e-8, int miter = 50,
                            bool verbose = false);

  static std::string type();
  static ParameterSet parameters();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int damage(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  virtual int ddamage_dd(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  virtual int ddamage_de(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
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

/// The standard damage model where the damage rate goes as the plastic strain
class NEMLStandardScalarDamagedModel_sd: public NEMLScalarDamagedModel_sd {
 public:
  NEMLStandardScalarDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::shared_ptr<NEMLModel_sd> base,
      std::shared_ptr<Interpolate> alpha = nullptr,
      double tol = 1.0e-8, int miter = 50,
      bool verbose = false);
  NEMLStandardScalarDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::shared_ptr<NEMLModel_sd> base,
      double alpha = 0.0,
      double tol = 1.0e-8, int miter = 50,
      bool verbose = false);

  virtual int damage(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  virtual int ddamage_dd(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  virtual int ddamage_de(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;
  virtual int ddamage_ds(double d_np1, double d_n, 
                     const double * const e_np1, const double * const e_n,
                     const double * const s_np1, const double * const s_n,
                     double T_np1, double T_n,
                     double t_np1, double t_n,
                     double * const dd) const;

  virtual int f(const double * const s_np1, double d_np1, 
                double T_np1, double & f) const = 0;
  virtual int df_ds(const double * const s_np1, double d_np1,
                   double T_np1, double * const df) const = 0;
  virtual int df_dd(const double * const s_np1, double d_np1,
                   double T_np1, double & df) const = 0;

 protected:
  double dep(const double * const s_np1, const double * const s_n,
             const double * const e_np1, const double * const e_n,
             double T_np1) const;

};

/// Simple power law damage
class NEMLPowerLawDamagedModel_sd: public NEMLStandardScalarDamagedModel_sd {
 public:
  NEMLPowerLawDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::shared_ptr<Interpolate> A, std::shared_ptr<Interpolate> a, 
      std::shared_ptr<NEMLModel_sd> base,
      std::shared_ptr<Interpolate> alpha = nullptr,
      double tol = 1.0e-8, int miter = 50,
      bool verbose = false);
  NEMLPowerLawDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      double A, double a,
      std::shared_ptr<NEMLModel_sd> base,
      double alpha = 0.0,
      double tol = 1.0e-8, int miter = 50,
      bool verbose = false);

  static std::string type();
  static ParameterSet parameters();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int f(const double * const s_np1, double d_np1,
                double T_np1, double & f) const;
  virtual int df_ds(const double * const s_np1, double d_np1, double T_np1,
                double * const df) const;
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
class NEMLExponentialWorkDamagedModel_sd: public NEMLStandardScalarDamagedModel_sd {
 public:
  NEMLExponentialWorkDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      std::shared_ptr<Interpolate> W0, std::shared_ptr<Interpolate> k0,
      std::shared_ptr<Interpolate> af,
      std::shared_ptr<NEMLModel_sd> base,
      std::shared_ptr<Interpolate> alpha = nullptr,
      double tol = 1.0e-8, int miter = 50,
      bool verbose = false);
  NEMLExponentialWorkDamagedModel_sd(
      std::shared_ptr<LinearElasticModel> elastic,
      double W0, double k0, double af,
      std::shared_ptr<NEMLModel_sd> base,
      double alpha = 0.0,
      double tol = 1.0e-8, int miter = 50,
      bool verbose = false);

  static std::string type();
  static ParameterSet parameters();
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  virtual int f(const double * const s_np1, double d_np1,
                double T_np1, double & f) const;
  virtual int df_ds(const double * const s_np1, double d_np1, double T_np1,
                double * const df) const;
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

} //namespace neml

#endif // DAMAGE_H
