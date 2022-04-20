#ifndef LARSONMILLER_H
#define LARSONMILLER_H

#include "interpolate.h"
#include "objects.h"
#include "solvers.h"

namespace neml {

class LMTrialState: public TrialState {
 public:
  virtual ~LMTrialState() {};
  double stress;
};

/// Classical Larson-Miller relation between stress and rupture time
//    Input function gives log(s) = fn(LMP) with LMP = T * (C + log(tr)) 
class NEML_EXPORT LarsonMillerRelation: public NEMLObject, public Solvable {
 public:
  LarsonMillerRelation(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  
  /// Stress as a function of rupture time (not really used)
  void sR(double t, double T, double & s) const;

  /// Rupture time as a function of stress
  void tR(double s, double T, double & t);

  /// Derivative of rupture time with respect to stress
  void dtR_ds(double s, double T, double & dt);

  /// Number of solver parameters
  virtual size_t nparams() const;
  /// Setup an iteration vector in the solver
  virtual void init_x(double * const x, TrialState * ts);
  /// Solver function returning the residual and jacobian of the nonlinear
  /// system of equations integrating the model
  virtual void RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J);

 protected:
  std::shared_ptr<Interpolate> fn_;
  double C_;
  double rtol_, atol_;
  int miter_;
  bool verbose_, linesearch_;
};

static Register<LarsonMillerRelation> regLarsonMillerRelation;

}

#endif // LARSONMILLER_H
