#ifndef SOLVERS_H
#define SOLVERS_H

#include <cstddef>
#include <memory>

#ifdef SOLVER_NOX
#include "NOX.H"
#include "NOX_LAPACK_Group.H"
#endif

namespace neml {

/// Generic nonlinear solver interface
class Solvable {
 public:
  virtual size_t nparams() const = 0;
  virtual int init_x(double * const x) = 0;
  virtual int RJ(const double * const x, double * const R, double * const J) = 0;
};

// This is entirely for testing
class TestRosenbrock: public Solvable {
 public:
  TestRosenbrock(size_t N);

  virtual size_t nparams() const;
  virtual int init_x(double * const x);
  virtual int RJ(const double * const x, double * const R, double * const J);

 private:
  const size_t N_;
};

/// Call the built-in solver
int solve(Solvable * system, double * x, 
          double tol = 1.0e-8, int miter = 50,
          bool verbose = false);

/// Default solver: plain NR
int newton(Solvable * system, double * x, 
          double tol, int miter, bool verbose);

#ifdef SOLVER_NOX
/// NOX OO interface
class NOXSolver: public NOX::LAPACK::Interface {
 public:
  NOXSolver(Solvable * system);

  const NOX::LAPACK::Vector& getInitialGuess();

  bool computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector& x);
  bool computeJacobian(NOX::LAPACK::Matrix<double>& J, 
      const NOX::LAPACK::Vector& x);

 private:
  NOX::LAPACK::Vector nox_guess_;
  std::shared_ptr<Solvable> system_;
};

/// Interface to nox
int nox(std::shared_ptr<Solvable> system, double * x, 
        double tol, int miter, bool verbose);

#endif

/// Helper to get numerical jacobian
int diff_jac(Solvable * system, const double * const x,
             double * const nJ, double eps = 1.0e-9);
/// Helper to get checksum
double diff_jac_check(Solvable * system, const double * const x,
                      const double * const J);

} // namespace neml

#endif // SOLVERS_H
