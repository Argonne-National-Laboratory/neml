#ifndef SINGLECRYSTAL_H
#define SINGLECRYSTAL_H

#include "kinematics.h"
#include "crystallography.h"

#include "../models.h"
#include "../elasticity.h"
#include "../history.h"

#include "../math/rotations.h"
#include "../math/tensors.h"

namespace neml {

/// Store the trial state for the single crystal model
class SCTrialState: public TrialState {
 public:

};

/// Single crystal model integrator
class SingleCrystalModel: public NEMLModel_ldi, public Solvable
{
 public:
  SingleCrystalModel(std::shared_ptr<KinematicModel> kinematics, 
                     std::shared_ptr<Orientation> initial_angle,
                     std::shared_ptr<Interpolate> alpha);
  virtual ~SingleCrystalModel();

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

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

   virtual size_t nhist() const;
   virtual int init_hist(double * const hist) const;

   virtual double alpha(double T) const;
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

 private:
  std::shared_ptr<KinematicModel> kinematics_;
  std::shared_ptr<Orientation> q0_;
  std::shared_ptr<Interpolate> alpha_;
};

static Register<SingleCrystalModel> regSingleCrystalModel;

} // namespace neml

#endif // SINGLECRYSTAL_H
