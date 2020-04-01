#ifndef SOLVERS_H
#define SOLVERS_H

#include <cstddef>
#include <memory>

#include "windows.h"

#ifdef SOLVER_NOX
#include "NOX.H"
#include "NOX_LAPACK_Group.H"
#endif

namespace neml {

/// Trial state
///  Store data the solver needs and can be passed into solution interface
class TrialState {
 public:
  virtual ~TrialState() {};
};

/// Generic nonlinear solver interface
class NEML_EXPORT Solvable {
 public:
  virtual ~Solvable() {};

  /// Number of parameters in the nonlinear equation
  virtual size_t nparams() const = 0;
  /// Initialize a guess to start the solution iterations
  virtual int init_x(double * const x, TrialState * ts) = 0;
  /// Nonlinear residual equations and corresponding jacobian
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J) = 0;
};

/// Call the built-in solver
int NEML_EXPORT solve(Solvable * system, double * x, TrialState * ts,
          double tol = 1.0e-8, int miter = 50,
          bool verbose = false, bool relative = false,
          double * R = nullptr, double * J = nullptr);

/// Default solver: plain NR
int NEML_EXPORT newton(Solvable * system, double * x, TrialState * ts,
          double tol, int miter, bool verbose, bool relative,
          double * R, double * J);

#ifdef SOLVER_NOX
/// NOX object-oriented interface
class NEML_EXPORT NOXSolver: public NOX::LAPACK::Interface {
 public:
  /// Setup with the solvable object and the trial state
  NOXSolver(Solvable * system, TrialState * ts);

  /// Get a NOX initial guess
  const NOX::LAPACK::Vector& getInitialGuess();

  /// Compute the residual
  bool computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector& x);
  /// Compute the jacobian
  bool computeJacobian(NOX::LAPACK::Matrix<double>& J,
      const NOX::LAPACK::Vector& x);

 private:
  NOX::LAPACK::Vector nox_guess_;
  Solvable * system_;
  TrialState * ts_;
};

/// Interface to nox
int NEML_EXPORT nox(Solvable * system, double * x, TrialState * ts,
        double tol, int miter, bool verbose, double * R, double * J);

#endif

/// Helper to get numerical jacobian
int NEML_EXPORT diff_jac(Solvable * system, const double * const x, TrialState * ts,
             double * const nJ, double eps = 1.0e-9);
/// Helper to get checksum
double NEML_EXPORT diff_jac_check(Solvable * system, const double * const x, TrialState * ts,
                      const double * const J);

} // namespace neml

#endif // SOLVERS_H
