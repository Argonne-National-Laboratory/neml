#include "solvers.h"

#include "nemlmath.h"
#include "nemlerror.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

namespace neml {

// This function is configured by the build
int solve(Solvable * system, double * x, TrialState * ts,
          double tol, int miter, bool verbose, bool relative)
{
#ifdef SOLVER_NOX
  return nox(system, x, ts, tol, miter, verbose);
#elif SOLVER_NEWTON
  // Actually selected the newton solver
  return newton(system, x, ts, tol, miter, verbose, relative);
#else
  // Default solver: plain NR
  return newton(system, x, ts, tol, miter, verbose, relative);
#endif
}

int newton(Solvable * system, double * x, TrialState * ts,
          double tol, int miter, bool verbose, bool relative)
{
  int n = system->nparams();
  system->init_x(x, ts);

  std::vector<double> Rv(n);
  std::vector<double> Jv(n*n);

  double * R = &Rv[0];
  double * J = &Jv[0];

  int ier = 0;

  ier = system->RJ(x, ts, R, J);
  if (ier != SUCCESS) return ier;

  double nR = norm2_vec(R, n);
  double nR0 = nR;
  int i = 0;

  if (verbose) {
    std::cout << "Iter.\tnR\t\tJe\t\tcn" << std::endl;
    double Jf = diff_jac_check(system, x, ts, J);
    double cn = condition(J, system->nparams());
    std::cout << std::setw(6) << std::left << i 
        << "\t" << std::setw(8) << std::left << std::scientific << nR 
        << "\t" << std::setw(8) << std::left << std::scientific << Jf
        << "\t" << std::setw(8) << std::left << std::scientific << cn
        << std::endl;
  }

  while ((nR > tol) && (i < miter))
  {
    if (relative) {
      if ((nR / nR0) < tol) break;
    }
    solve_mat(J, n, R);

    for (int j=0; j<n; j++) x[j] -= R[j];

    system->RJ(x, ts, R, J);
    nR = norm2_vec(R, n);
    i++;

    if (verbose) {
      double Jf = diff_jac_check(system, x, ts, J);
      double cn = condition(J, system->nparams());
      std::cout << i << "\t" << nR << "\t" << Jf << "\t" << cn << std::endl;
    }
  }

  if (verbose) {
    std::cout << std::endl;
  }

  if (ier != SUCCESS) return ier;

  if (i == miter) return MAX_ITERATIONS;

  return SUCCESS;
}

/// Helper to get numerical jacobian
int diff_jac(Solvable * system, const double * const x, TrialState * ts,
             double * const nJ, double eps)
{
  std::vector<double> R0v(system->nparams());
  std::vector<double> nRv(system->nparams());
  std::vector<double> nXv(system->nparams());
  std::vector<double> dJv(system->nparams() * system->nparams());

  double * R0 = &R0v[0];
  double * nR = &nRv[0];
  double * nX = &nXv[0];
  double * dJ = &dJv[0];

  system->RJ(x, ts, R0, dJ);
  
  for (size_t i=0; i<system->nparams(); i++) {
    std::copy(x, x+system->nparams(), nX);
    double dx = eps * fabs(nX[i]);
    if (dx < eps) dx = eps;
    nX[i] += dx;
    system->RJ(nX, ts, nR, dJ);
    for (size_t j=0; j<system->nparams(); j++) {
      nJ[CINDEX(j,i,system->nparams())] = (nR[j] - R0[j]) / dx;
    }
  }
  
  return 0;
}

/// Helper to get checksum
double diff_jac_check(Solvable * system, const double * const x,
                      TrialState * ts, const double * const J)
{
  std::vector<double> nJv(system->nparams() * system->nparams());
  double * nJ = &nJv[0];
  
  diff_jac(system, x, ts, nJ);
  double ss = 0.0;
  double js = 0.0;
  for (size_t i=0; i< system->nparams() * system->nparams(); i++) {
    ss += pow(J[i] - nJ[i], 2.0);
    js += pow(J[i], 2.0);
  }

  return ss/js;
}

// START NOX STUFF
#ifdef SOLVER_NOX
NOXSolver::NOXSolver(Solvable * system) :
    system_(system), nox_guess_(system->nparams()), ts_(ts)
{
  std::vector<double> xn(system_->nparams());
  double * x = &xn[0];
  system_->init_x(x, ts_);
  for (int i=0; i<system_->nparams(); i++) {
    nox_guess_(i) = x[i];
  }
}

const NOX::LAPACK::Vector& NOXSolver::getInitialGuess()
{
  return nox_guess_;
}

bool NOXSolver::computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector& x)
{
  // This is highly inefficient
  std::vector<double> Riv(system_->nparams());
  std::vector<double> Jiv(system_->nparams());
  std::vector<double> xiv(system_->nparams());
  
  double * Ri = &Riv[0];
  double * Ji = &Jiv[0];
  double * xi = &xiv[0];

  for (int i=0; i<system_->nparams(); i++) {
    xi[i] = x(i);
  }
  system_->RJ(xi, ts_, Ri, Ji);
  
  for (int i=0; i<system_->nparams(); i++) {
    f(i) = Ri[i];
  }

  return true;
}

bool NOXSolver::computeJacobian(NOX::LAPACK::Matrix<double>& J,
                                const NOX::LAPACK::Vector & x)
{
  // This is highly inefficient
  std::vector<double> Riv(system_->nparams());
  std::vector<double> Jiv(system_->nparams());
  std::vector<double> xiv(system_->nparams());
  
  double * Ri = &Riv[0];
  double * Ji = &Jiv[0];
  double * xi = &xiv[0];

  for (int i=0; i<system_->nparams(); i++) {
    xi[i] = x(i);
  }
  system_->RJ(xi, ts_, Ri, Ji);

  for (int i=0; i<system_->nparams(); i++) {
    for (int j=0; j<system_->nparams(); j++) {
      J(i,j) = Ji[CINDEX(i,j,system_->nparams())];
    }
  }

  return true;
}


int nox(Solvable * system, double * x, TrialState * ts,
        double tol, int miter, bool verbose)
{
  // Setup solver
  NOXSolver solver(system, ts);

  // Setup NOX
  Teuchos::RCP<NOX::LAPACK::Group> grp = 
      Teuchos::rcp(new NOX::LAPACK::Group(solver));

  // Setup status tests
  Teuchos::RCP<NOX::StatusTest::NormF> statusTestF = 
      Teuchos::rcp(new NOX::StatusTest::NormF(tol));
  Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestI = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(miter));
  Teuchos::RCP<NOX::StatusTest::Combo> status = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                              statusTestF, statusTestI));

  // Setup solver parameters
  Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr = 
      Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& solverParameters = *solverParametersPtr;

  // Verbosity
  if (verbose) {
    solverParameters.sublist("Printing").set("Output Information",
                                             NOX::Utils::Warning +
                                             NOX::Utils::OuterIteration);
  }
  else {
    solverParameters.sublist("Printing").set("Output Information",
                                             0);
  }

  // Actual solver
  solverParameters.set("Nonlinear Solver", "Line Search Based");
  Teuchos::ParameterList& lineSearchParameters = 
      solverParameters.sublist("Line Search");
  lineSearchParameters.set("Method", "Backtrack");
  

  Teuchos::ParameterList& directionParameters = 
      solverParameters.sublist("Direction");
  directionParameters.set("Method", "Newton");

  // Actual solver
  Teuchos::RCP<NOX::Solver::Generic> nox_solver = 
      NOX::Solver::buildSolver(grp, status, solverParametersPtr);

  // Actual solve
  NOX::StatusTest::StatusType result = nox_solver->solve();

  // Check if we actually succeeded
  if (result != NOX::StatusTest::Converged) {
    return MAX_ITERATIONS; 
  }

  // Get the solution
  NOX::LAPACK::Group solnGrp = 
      dynamic_cast<const NOX::LAPACK::Group&>(nox_solver->getSolutionGroup());
  NOX::LAPACK::Vector soln = dynamic_cast<const NOX::LAPACK::Vector&>(
      solnGrp.getX());
  for (int i=0; i<system->nparams(); i++) {
    x[i] = soln(i);
  }

  return 0;
}

#endif
} // namespace neml
