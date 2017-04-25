#ifndef CREEP_H
#define CREEP_H

#include "interpolate.h"
#include "solvers.h"

namespace neml {

/// Scalar creep functions in terms of effective stress and strain
class ScalarCreepRule {
  public:
   virtual int g(double seq, double eeq, double t, double T, double & g) const = 0;
   virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const = 0;
   virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const = 0;
   virtual int dg_dt(double seq, double eeq, double t, double T, double & dg) const;
   virtual int dg_dT(double seq, double eeq, double t, double T, double & dg) const;
};

/// Simple power law creep
class PowerLawCreep: public ScalarCreepRule {
  public:
   PowerLawCreep(double A, double n);
   PowerLawCreep(std::shared_ptr<Interpolate> A, std::shared_ptr<Interpolate> n);

   virtual int g(double seq, double eeq, double t, double T, double & g) const;
   virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const;
   virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const;

   double A(double T) const;
   double n(double T) const;

  private:
   const std::shared_ptr<const Interpolate> A_, n_;
};

/// Classical Norton-Bailey creep
class NortonBaileyCreep: public ScalarCreepRule {
 public:
  NortonBaileyCreep(double A, double m, double n);
  NortonBaileyCreep(std::shared_ptr<Interpolate> A, std::shared_ptr<Interpolate> m,
                    std::shared_ptr<Interpolate> n);

   virtual int g(double seq, double eeq, double t, double T, double & g) const;
   virtual int dg_ds(double seq, double eeq, double t, double T, double & dg) const;
   virtual int dg_de(double seq, double eeq, double t, double T, double & dg) const;

   double A(double T) const;
   double m(double T) const;
   double n(double T) const;

 private:
   const std::shared_ptr<const Interpolate> A_, m_, n_;

};

/// Creep trial state
class CreepModelTrialState : public TrialState {
 public:
  double T, dt, t;
  double s_np1[6];
  double e_n[6];
};

/// Master class of all creep models defining the interface
class CreepModel: public Solvable {
 public:
  int update(const double * const s_np1, 
             double * const e_np1, const double * const e_n,
             double T_np1, double T_n,
             double t_np1, double t_n,
             double * const A_np1) const;

  virtual int f(const double * const s, const double * const e, double t, double T, 
                double * const f) const = 0;
  virtual int df_ds(const double * const s, const double * const e, double t, double T, 
                double * const df) const = 0;
  virtual int df_de(const double * const s, const double * const e, double t, double T, 
                double * const df) const = 0;
  virtual int df_dt(const double * const s, const double * const e, double t, double T, 
                double * const df) const ;
  virtual int df_dT(const double * const s, const double * const e, double t, double T, 
                double * const df) const;
  
  // Solvable implementation
  int make_trial_state(const double * const s_np1, 
                       const double * const e_n,
                       double T_np1, double T_n,
                       double t_np1, double t_n,
                       CreepModelTrialState & ts) const;
  virtual size_t nparams() const;
  virtual int init_x(double * const x, TrialState * ts);
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

};

/// J2 creep based on a scalar creep rule
class J2CreepModel: public CreepModel {
 public:
  J2CreepModel(std::shared_ptr<ScalarCreepRule> rule);

  virtual int f(const double * const s, const double * const e, double t, double T, 
                double * const f) const;
  virtual int df_ds(const double * const s, const double * const e, double t, double T, 
                double * const df) const;
  virtual int df_de(const double * const s, const double * const e, double t, double T, 
                double * const df) const;
  virtual int df_dt(const double * const s, const double * const e, double t, double T, 
                double * const df) const ;
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

} // namespace neml

#endif // CREEP_H
