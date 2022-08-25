#ifndef SINGLECRYSTAL_H
#define SINGLECRYSTAL_H

#include "kinematics.h"
#include "crystallography.h"
#include "postprocessors.h"

#include "../models.h"
#include "../elasticity.h"
#include "../history.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

#include "../windows.h"

namespace neml {

class CrystalPostprocessor; // forward declaration

/// Store the trial state for the single crystal model
class SCTrialState: public TrialState {
 public:
  SCTrialState(const Symmetric & d, const Skew & w, const Symmetric & S, 
               const Symmetric & S_n, const History & H, 
               const Orientation & Q, const Lattice & lattice,
               double T, double dt,
               const History & fixed) :
      d(d), w(w), S(S), S_n(S_n), history(H), Q(Q), lattice(lattice), T(T), dt(dt),
      fixed(fixed)
  {};

  Symmetric d;
  Skew w;
  Symmetric S;
  Symmetric S_n;
  History history;
  Orientation Q;
  Lattice lattice;
  double T;
  double dt;
  History fixed;
};

/// Single crystal model integrator
class NEML_EXPORT SingleCrystalModel: public NEMLModel_ldi, public Solvable
{
 public:
  /// Raw constructor
  SingleCrystalModel(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Setup blank history
  void populate_state(History & history) const;
  /// Actually initialize history
  void init_state(History & history) const;
  
  /// Not changing state
  void populate_static(History & history) const;
  /// Static stat
  void init_static(History & history) const;

  /// Useful methods for external models that want an idea of an average
  /// strength
  double strength(double * const hist, double T) const;

  /// Used to calculate the Nye tensor
  void Fe(double * const stress, double * const hist, double T,
          double * Fe) const;
  
  /// Large deformation incremental update
  virtual void update_ld_inc(
       const double * const d_np1, const double * const d_n,
       const double * const w_np1, const double * const w_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1, double * const B_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n);

  /// Instantaneous CTE
  virtual double alpha(double T) const;

  /// Helper to calculate the elastic strain
  virtual void elastic_strains(const double * const s_np1,
                               double T_np1, const double * const h_np1,
                               double * const e_np1) const;

  /// Number of nonlinear equations to solve in the integration
  virtual size_t nparams() const;
  /// Setup an initial guess for the nonlinear solution
  virtual void init_x(double * const x, TrialState * ts);
  /// Integration residual and jacobian equations
  virtual void RJ(const double * const x, TrialState * ts, double * const R,
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
  void set_active_orientation(double * const hist, const Orientation & q);
  /// Set the current orientation given an active rotation (crystal to lab)
  void set_active_orientation(History & hist, const Orientation & q);
  /// Set the current orientation given a passive rotation (lab to crystal)
  void set_passive_orientation(double * const hist, const Orientation & q);
  /// Set the current orientation given a passive rotation (lab to crystal)
  void set_passive_orientation(History & hist, const Orientation & q);

  /// Whether this model uses the nye tensor
  virtual bool use_nye() const;

  /// Actually update the Nye tensor
  void update_nye(double * const hist, const double * const nye) const;

 private:
  void attempt_update_ld_inc_(
       const double * const d_np1, const double * const d_n,
       const double * const w_np1, const double * const w_n,
       double T_np1, double T_n,
       double t_np1, double t_n,
       double * const s_np1, const double * const s_n,
       double * const h_np1, const double * const h_n,
       double * const A_np1, double * const B_np1,
       double & u_np1, double u_n,
       double & p_np1, double p_n, int trial_type);

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

  void solve_substep_(SCTrialState * ts, Symmetric & stress, History & hist);

 private:
  std::shared_ptr<KinematicModel> kinematics_;
  std::shared_ptr<Lattice> lattice_;
  std::shared_ptr<Orientation> q0_;
  std::shared_ptr<Interpolate> alpha_;
  bool update_rotation_;
  double rtol_, atol_;
  int miter_;
  bool verbose_, linesearch_;
  int max_divide_;

  std::vector<std::shared_ptr<CrystalPostprocessor>> postprocessors_;

  bool elastic_predictor_, fallback_elastic_predictor_;
  int force_divide_;
  bool elastic_predictor_first_step_;
};

static Register<SingleCrystalModel> regSingleCrystalModel;

} // namespace neml

#endif // SINGLECRYSTAL_H
