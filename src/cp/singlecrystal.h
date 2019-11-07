#ifndef SINGLECRYSTAL_H
#define SINGLECRYSTAL_H

#include "kinematics.h"
#include "crystallography.h"

#include "../models.h"
#include "../elasticity.h"
#include "../history.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include "../windows.h"

namespace neml {

/// Store the trial state for the single crystal model
class SCTrialState: public TrialState {
 public:
  SCTrialState(const Symmetric & d, const Skew & w, const Symmetric & S, const
               History & H, const Orientation & Q, const Lattice & lattice,
               double T, double dt,
               const Symmetric & s_guess, const History & h_guess,
               const History & fixed) :
      d(d), w(w), S(S), history(H), Q(Q), lattice(lattice), T(T), dt(dt),
      s_guess(s_guess), h_guess(h_guess), fixed(fixed)
  {};

  Symmetric d;
  Skew w;
  Symmetric S;
  History history;
  Orientation Q;
  Lattice lattice;
  double T;
  double dt;
  Symmetric s_guess;
  History h_guess;
  History fixed;
};

/// Single crystal model integrator
class SingleCrystalModel: public NEMLModel_ldi, public Solvable
{
 public:
  /// Raw constructor
  SingleCrystalModel(std::shared_ptr<KinematicModel> kinematics,
                     std::shared_ptr<Lattice> lattice,
                     std::shared_ptr<Orientation> initial_angle,
                     std::shared_ptr<Interpolate> alpha,
                     bool update_rotation,
                     double tol, int miter, bool verbose,
                     int max_divide);
  /// Destructor
  virtual ~SingleCrystalModel();

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Setup blank history
  void populate_history(History & history) const;

  /// Actually initialize history
  void init_history(History & history) const;

  /// Useful methods for external models that want an idea of an average
  /// strength
  double strength(const History & history, Lattice & L, double T) const;

  /// Large deformation incremental update
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

  /// Number of stored history variables
  virtual size_t nhist() const;
  /// Initialize history raw pointer array
  virtual int init_hist(double * const hist) const;

  /// Instantaneous CTE
  virtual double alpha(double T) const;

  /// Helper to calculate the elastic strain
  virtual int elastic_strains(const double * const s_np1,
                               double T_np1, const double * const h_np1,
                               double * const e_np1) const;

  /// Number of nonlinear equations to solve in the integration
  virtual size_t nparams() const;
  /// Setup an initial guess for the nonlinear solution
  virtual int init_x(double * const x, TrialState * ts);
  /// Integration residual and jacobian equations
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

  /// Get the current orientation in the active convention (raw ptr history)
  Orientation get_active_orientation(double * const hist) const;
  /// Get the current orientation in the active convention
  Orientation get_active_orientation(const History & hist) const;
  /// Get the current orientation in the passive convention (raw ptr history)
  Orientation get_passive_orientation(double * const hist) const;
  /// Get the current orientation in the passive convention
  Orientation get_passive_orientation(const History & hist) const;

  /// Set the current orientation given an active rotation (crystal to lab)
  void set_active_orientation(double * const hist, const Orientation & q) const;
  /// Set the current orientation given an active rotation (crystal to lab)
  void set_active_orientation(History & hist, const Orientation & q) const;
  /// Set the current orientation given a passive rotation (lab to crystal)
  void set_passive_orientation(double * const hist, const Orientation & q) const;
  /// Set the current orientation given a passive rotation (lab to crystal)
  void set_passive_orientation(History & hist, const Orientation & q) const;

 private:
  History gather_history_(double * data) const;
  History gather_history_(const double * data) const;
  History gather_blank_history_() const;

  void calc_tangents_(Symmetric & S, History & H,
                      SCTrialState * ts, double * const A, double * const B);
  Orientation update_rot_(Symmetric & S_np1, History & H_np1, SCTrialState * ts) const;
  double calc_energy_inc_(const Symmetric & D_np1, const Symmetric & D_n,
                          const Symmetric & S_np1, const Symmetric & S_n) const;
  double calc_work_inc_(const Symmetric & D_np1, const Symmetric & D_n,
                        const Symmetric & S_np1, const Symmetric & S_n,
                        double T_np1, double T_n, const Orientation & Q_np1,
                        const Orientation & Q_n, const History & H_np1,
                        const History & H_n) const;

  int solve_substep_(SCTrialState * ts, Symmetric & stress, History & hist);

 private:
  std::shared_ptr<KinematicModel> kinematics_;
  std::shared_ptr<Lattice> lattice_;
  std::shared_ptr<Orientation> q0_;
  std::shared_ptr<Interpolate> alpha_;
  bool update_rotation_;
  double tol_;
  int miter_;
  bool verbose_;
  int max_divide_;

  History stored_hist_;
};

static Register<SingleCrystalModel> regSingleCrystalModel;

} // namespace neml

#endif // SINGLECRYSTAL_H
