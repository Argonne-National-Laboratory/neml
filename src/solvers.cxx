#include "solvers.h"

#include "nemlmath.h"
#include "nemlerror.h"

#include <algorithm>
#include <iostream>
#include <iomanip>

namespace neml {

TestRosenbrock::TestRosenbrock(size_t N) :
    N_(N)
{

}

size_t TestRosenbrock::nparams() const
{
  return N_;
}

int TestRosenbrock::init_x(double * const x, TrialState * ts)
{
  std::fill(x, x+N_, 0.25);
  return 0;
}

int TestRosenbrock::RJ(const double * const x, TrialState * ts, 
                       double * const R, double * const J)
{
  R[0] = 1.0 - x[0];
  for (int i=1; i<N_; i++) {
    R[i] = 10.0 * (x[i] - x[i-1] * x[i-1]);
  }

  std::fill(J, J+(N_*N_), 0.0);
  J[0] = -1.0;
  for (int i=1; i<N_; i++) {
    J[CINDEX(i,i-1,N_)] = -20.0 * x[i-1];
    J[CINDEX(i,i,N_)] = 10.0;
  }

  return 0;
}


// This function is configured by the build
int solve(Solvable * system, double * x, TrialState * ts,
          double tol, int miter, bool verbose)
{
#ifdef SOLVER_NOX
  return nox(system, x, ts, tol, miter, verbose);
#elif SOLVER_NEWTON
  // Actually selected the newton solver
  return newton(system, x, ts, tol, miter, verbose);
#else
  // Default solver: plain NR
  return newton(system, x, ts, tol, miter, verbose);
#endif
}

int newton(Solvable * system, double * x, TrialState * ts,
          double tol, int miter, bool verbose)
{
  int n = system->nparams();
  system->init_x(x, ts);
  double * R = new double[n];
  double * J = new double[n*n];
  
  system->RJ(x, ts, R, J);

  double nR = norm2_vec(R, n);
  int i = 0;
  int ier = 0;

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

  delete [] R;
  delete [] J;

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
  double * R0 = new double[system->nparams()];
  double * nR = new double[system->nparams()];
  double * nX = new double[system->nparams()];
  double * dJ = new double[system->nparams()*system->nparams()];

  system->RJ(x, ts, R0, dJ);
  
  for (int i=0; i<system->nparams(); i++) {
    std::copy(x, x+system->nparams(), nX);
    double dx = eps * fabs(nX[i]);
    if (dx < eps) dx = eps;
    nX[i] += dx;
    system->RJ(nX, ts, nR, dJ);
    for (int j=0; j<system->nparams(); j++) {
      nJ[CINDEX(j,i,system->nparams())] = (nR[j] - R0[j]) / dx;
    }
  }
  
  delete [] R0;
  delete [] nR;
  delete [] nX;
  delete [] dJ;
  return 0;
}

/// Helper to get checksum
double diff_jac_check(Solvable * system, const double * const x,
                      TrialState * ts, const double * const J)
{
  double * nJ = new double[system->nparams() * system->nparams()];
  
  diff_jac(system, x, ts, nJ);
  double ss = 0.0;
  double js = 0.0;
  for (int i=0; i< system->nparams() * system->nparams(); i++) {
    ss += pow(J[i] - nJ[i], 2.0);
    js += pow(J[i], 2.0);
  }

  delete [] nJ;

  return ss/js;
}

// START NOX STUFF
#ifdef SOLVER_NOX
NOXSolver::NOXSolver(Solvable * system) :
    system_(system), nox_guess_(system->nparams()), ts_(ts)
{
  double * x = new double[system_->nparams()];
  system_->init_x(x, ts_);
  for (int i=0; i<system_->nparams(); i++) {
    nox_guess_(i) = x[i];
  }
  delete [] x;
}

const NOX::LAPACK::Vector& NOXSolver::getInitialGuess()
{
  return nox_guess_;
}

bool NOXSolver::computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector& x)
{
  // This is highly inefficient
  double * Ri = new double [system_->nparams()];
  double * Ji = new double [system_->nparams() * system_->nparams()];
  double * xi = new double [system_->nparams()];
  for (int i=0; i<system_->nparams(); i++) {
    xi[i] = x(i);
  }
  system_->RJ(xi, ts_, Ri, Ji);
  
  for (int i=0; i<system_->nparams(); i++) {
    f(i) = Ri[i];
  }

  delete [] Ri;
  delete [] Ji;
  delete [] xi;

  return true;
}

bool NOXSolver::computeJacobian(NOX::LAPACK::Matrix<double>& J,
                                const NOX::LAPACK::Vector & x)
{
  // This is highly inefficient
  double * Ri = new double [system_->nparams()];
  double * Ji = new double [system_->nparams() * system_->nparams()];
  double * xi = new double [system_->nparams()];
  for (int i=0; i<system_->nparams(); i++) {
    xi[i] = x(i);
  }
  system_->RJ(xi, ts_, Ri, Ji);

  for (int i=0; i<system_->nparams(); i++) {
    for (int j=0; j<system_->nparams(); j++) {
      J(i,j) = Ji[CINDEX(i,j,system_->nparams())];
    }
  }

  delete [] Ri;
  delete [] Ji;
  delete [] xi;

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
  
  //solverParameters.set("Nonlinear Solver", "Trust Region Based");
  //Teuchos::ParameterList& lineSearchParameters = 
  //    solverParameters.sublist("Trust Region");


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
