#ifndef SOLVERS_H
#define SOLVERS_H

#include <cstddef>
#include <memory>

#ifdef SOLVER_NOX
#include "NOX.H"
#include "NOX_LAPACK_Group.H"
#endif

namespace neml {

/// Trial state classes
//  Store data the solver needs and can be passed into solution interface

class TrialState {

};

/// Generic nonlinear solver interface
class Solvable {
 public:
  virtual size_t nparams() const = 0;
  virtual int init_x(double * const x, TrialState * ts) = 0;
  virtual int RJ(const double * const x, TrialState * ts, double * const R,
                 double * const J) = 0;
};

/// Call the built-in solver
int solve(Solvable * system, double * x, TrialState * ts, 
          double tol = 1.0e-8, int miter = 50,
          bool verbose = false, bool relative = false);

/// Default solver: plain NR
int newton(Solvable * system, double * x, TrialState * ts,
          double tol, int miter, bool verbose, bool relative);

#ifdef SOLVER_NOX
/// NOX OO interface
class NOXSolver: public NOX::LAPACK::Interface {
 public:
  NOXSolver(Solvable * system, TrialState * ts);

  const NOX::LAPACK::Vector& getInitialGuess();

  bool computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector& x);
  bool computeJacobian(NOX::LAPACK::Matrix<double>& J, 
      const NOX::LAPACK::Vector& x);

 private:
  NOX::LAPACK::Vector nox_guess_;
  Solvable * system_;
  TrialState * ts_;
};

/// Interface to nox
int nox(Solvable * system, double * x, TrialState * ts, 
        double tol, int miter, bool verbose);

#endif

/// Helper to get numerical jacobian
int diff_jac(Solvable * system, const double * const x, TrialState * ts,
             double * const nJ, double eps = 1.0e-9);
/// Helper to get checksum
double diff_jac_check(Solvable * system, const double * const x, TrialState * ts,
                      const double * const J);

} // namespace neml

#endif // SOLVERS_H
