#ifndef DAMAGE_H
#define DAMAGE_H

#include "neml.h"
#include "elasticity.h"

namespace neml {

/// Small strain damage model
class NEMLDamagedModel_sd: public NEMLModel_sd {
 public:
  NEMLDamagedModel_sd(std::shared_ptr<NEMLModel_sd> base,
                      std::shared_ptr<Interpolate> alpha = nullptr);
  NEMLDamagedModel_sd(std::shared_ptr<NEMLModel_sd> base,
                      double alpha = 0.0);
  virtual ~NEMLDamagedModel_sd();

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
  virtual int elastic_strains(const double * const s_np1,
                              double T_np1,
                              double * const e_np1) const = 0;

  virtual size_t ndamage() const = 0;
  virtual int init_damage(double * const damage) const = 0;

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
  NEMLScalarDamagedModel_sd(std::shared_ptr<NEMLModel_sd> base,
                            std::shared_ptr<Interpolate> alpha = nullptr);
  NEMLScalarDamagedModel_sd(std::shared_ptr<NEMLModel_sd> base,
                            double alpha = 0.0);
  virtual ~NEMLScalarDamagedModel_sd();

  virtual int update_sd(
      const double * const e_np1, const double * const e_n,
      double T_np1, double T_n,
      double t_np1, double t_n,
      double * const s_np1, const double * const s_n,
      double * const h_np1, const double * const h_n,
      double * const A_np1,
      double & u_np1, double u_n,
      double & p_np1, double p_n);
  virtual int elastic_strains(const double * const s_np1,
                              double T_np1,
                              double * const e_np1) const;

  virtual size_t ndamage() const;
  virtual int init_damage(double * const damage) const;

  virtual size_t nparams() const;
  virtual int init_x(double * const x, TrialState * ts);
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

  virtual int damage(double d, const double * const e_np1, 
                     double * const s_np1, double * const d_np1) const = 0;
  virtual int ddamage_dd(double d, const double * const e_np1, 
                     double * const s_np1, double * const dd) const = 0;
  virtual int ddamage_de(double d, const double * const e_np1, 
                     double * const s_np1, double * const dd) const = 0;
  virtual int ddamage_ds(double d, const double * const e_np1, 
                     double * const s_np1, double * const dd) const = 0;
};

/// The standard damage model where the damage rate goes as the plastic strain
class NEMLStandardDamagedModel_sd: public NEMLScalarDamagedModel_sd {


};

/// Simple power law damage
class NEMLPowerLawDamagedModel_sd: public NEMLStandardDamagedModel_sd {


};

} //namespace neml

#endif // DAMAGE_H
