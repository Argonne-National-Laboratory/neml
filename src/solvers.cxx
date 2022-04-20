#include "solvers.h"

#include "math/nemlmath.h"
#include "nemlerror.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

namespace neml {

// This function is configured by the build
void solve(Solvable * system, double * x, TrialState * ts,
          SolverParameters p, double * R, double * J)
{
#ifdef SOLVER_NOX
  nox(system, x, ts, p.atol, p.miter, p.verbose, R, J);
#elif SOLVER_NEWTON
  // Actually selected the newton solver
  newton(system, x, ts, p, R, J);
#else
  // Default solver: plain NR
  newton(system, x, ts, p, R, J);
#endif
}

void newton(Solvable * system, double * x, TrialState * ts, SolverParameters p, double * R,
           double * J)
{
  int mline = 10;

  int n = system->nparams();
  system->init_x(x, ts);

  bool local_R = (R==nullptr);
  bool local_J = (J==nullptr);

  if (local_R) {
    R = new double [n];
  }
  if (local_J) {
    J = new double [n*n];
  }

  system->RJ(x, ts, R, J);

  double nR = norm2_vec(R, n);
  double nR0 = nR;
  int i = 0;
  double alpha = 1.0;

  if (p.verbose) {
    std::cout << "Iter.\tnR\t\tJe\t\tcn" ;
    if (p.linesearch) std::cout << "\t\talpha";
    std::cout << std::endl;
    double Jf = diff_jac_check(system, x, ts, J);
    double cn = condition(J, system->nparams());
    std::cout << std::setw(6) << std::left << i 
        << "\t" << std::setw(8) << std::left << std::scientific << nR 
        << "\t" << std::setw(8) << std::left << std::scientific << Jf
        << "\t" << std::setw(8) << std::left << std::scientific << cn;
    if (p.linesearch) std::cout << "\t" << std::setw(8) << std::left << std::scientific << alpha;

    std::cout << std::endl;
  }

  while (true)
  {
    if ((nR < p.atol) || ((nR / nR0) < p.rtol)) break;

    solve_mat(J, n, R);

    if (p.linesearch) {
      int nsearch = 0;
      alpha = 1.0;
      double * x_orig = new double [n];
      std::copy(x, x+n, x_orig);
      double * dir = new double [n];
      std::copy(R, R+n, dir);
      double nRt = 0.0;
      bool linesearch_error = false;
      while (nsearch < mline) {
        for (int j=0; j<n; j++) x[j] = x_orig[j] - alpha * dir[j];
        system->RJ(x, ts, R, J);
        nRt = norm2_vec(R, n);
        if (nRt < nR) break;
        alpha /= 2.0;
        nsearch += 1;
      }
      delete [] x_orig;
      delete [] dir;
      if (linesearch_error) {
        break;
      }
      nR = nRt;
      if (nsearch == mline) {
        // Arguably it's better to let it fall through, but we could
        // raise a NonlinearSolverError right now with a line search
        // failed message
      }
    }
    else {
      for (int j=0; j<n; j++) x[j] -= R[j];
      system->RJ(x, ts, R, J);
      nR = norm2_vec(R, n);
    }
    i++;

    if (p.verbose) {
      double Jf = diff_jac_check(system, x, ts, J);
      double cn = condition(J, system->nparams());
      std::cout << i << "\t" << nR << "\t" << Jf << "\t" << cn;
      if (p.linesearch) std::cout << "\t" << alpha;
      std::cout << std::endl;
    }

    if (i >= p.miter) break;
  }

  if (p.verbose) {
    std::cout << std::endl;
  }

  if (local_R) {
    delete [] R;
  }

  if (local_J) {
    delete [] J;
  }

  if (i == p.miter) 
    throw NonlinearSolverError("Nonlinear solver exceeded maximum allowed iterations!");
}

/// Helper to get numerical jacobian
void diff_jac(Solvable * system, const double * const x, TrialState * ts,
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
NOXSolver::NOXSolver(Solvable * system, TrialState * ts) :
    nox_guess_(system->nparams()), system_(system), ts_(ts)
{
  std::vector<double> xn(system_->nparams());
  double * x = &xn[0];
  system_->init_x(x, ts_);
  for (size_t i=0; i<system_->nparams(); i++) {
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
  std::vector<double> Jiv(system_->nparams()*system_->nparams());
  std::vector<double> xiv(system_->nparams());
  
  double * Ri = &Riv[0];
  double * Ji = &Jiv[0];
  double * xi = &xiv[0];

  for (size_t i=0; i<system_->nparams(); i++) {
    xi[i] = x(i);
  }
  system_->RJ(xi, ts_, Ri, Ji);
  
  for (size_t i=0; i<system_->nparams(); i++) {
    f(i) = Ri[i];
  }

  return true;
}

bool NOXSolver::computeJacobian(NOX::LAPACK::Matrix<double>& J,
                                const NOX::LAPACK::Vector & x)
{
  // This is highly inefficient
  std::vector<double> Riv(system_->nparams());
  std::vector<double> Jiv(system_->nparams()*system_->nparams());
  std::vector<double> xiv(system_->nparams());
  
  double * Ri = &Riv[0];
  double * Ji = &Jiv[0];
  double * xi = &xiv[0];

  for (size_t i=0; i<system_->nparams(); i++) {
    xi[i] = x(i);
  }
  system_->RJ(xi, ts_, Ri, Ji);

  for (size_t i=0; i<system_->nparams(); i++) {
    for (size_t j=0; j<system_->nparams(); j++) {
      J(i,j) = Ji[CINDEX(i,j,system_->nparams())];
    }
  }

  return true;
}


void nox(Solvable * system, double * x, TrialState * ts,
        double tol, int miter, bool verbose, double * R,
        double * J)
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
    throw NonlinearSolverError("Nonlinear solver exceeded maximum allowed iterations!");
  }

  // Get the solution
  NOX::LAPACK::Group solnGrp = 
      dynamic_cast<const NOX::LAPACK::Group&>(nox_solver->getSolutionGroup());
  NOX::LAPACK::Vector soln = dynamic_cast<const NOX::LAPACK::Vector&>(
      solnGrp.getX());
  for (size_t i=0; i<system->nparams(); i++) {
    x[i] = soln(i);
  }
  
  // Would need to be figured out to be efficient
  if ((R != nullptr) || (J != nullptr)) {
    bool local_R = (R == nullptr);
    if (local_R) {
      R = new double [system->nparams()];
    }
    if (local_J) {
      J = new double [system->nparams() * system->nparams()];
    }

    system->RJ(x, ts, R, J);

    if (local_R) {
      delete [] R;
    }
    if (local_J) {
      delete [] J;
    }
  }
}

#endif

TestPower::TestPower(double A, double n, double b, double x0) :
    A_(A), n_(n), b_(b), x0_(x0)
{

}

size_t TestPower::nparams() const
{
  return 1;
}

void TestPower::init_x(double * const x, TrialState * ts)
{
  x[0] = x0_;
}

void TestPower::RJ(const double * const x, TrialState * ts, double * const R, 
       double * const J)
{
  R[0] = A_ * std::pow(x[0], n_) + b_;
  J[0] = A_ * n_ * std::pow(x[0], n_-1.0);
}


} // namespace neml
