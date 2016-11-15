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

int TestRosenbrock::init_x(double * const x)
{
  std::fill(x, x+N_, 0.25);
  return 0;
}

int TestRosenbrock::RJ(const double * const x, double * const R, 
                  double * const J)
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
int solve(std::shared_ptr<Solvable> system, double * x,
          double tol, int miter, bool verbose)
{
  // For now just plain NR
  return newton(system, x, tol, miter, verbose);
}

int newton(std::shared_ptr<Solvable> system, double * x,
          double tol, int miter, bool verbose)
{
  int n = system->nparams();
  system->init_x(x);
  double * R = new double[n];
  double * J = new double[n*n];
  
  system->RJ(x, R, J);

  double nR = norm2_vec(R, n);
  int i = 0;
  int ier = 0;

  if (verbose) {
    std::cout << "Iter.\tnR\t\tJe\t\tcn" << std::endl;
    double Jf = diff_jac_check(system, x, J);
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

    system->RJ(x, R, J);
    nR = norm2_vec(R, n);
    i++;

    if (verbose) {
      double Jf = diff_jac_check(system, x, J);
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
int diff_jac(std::shared_ptr<Solvable> system, const double * const x,
             double * const nJ, double eps)
{
  double * R0 = new double[system->nparams()];
  double * nR = new double[system->nparams()];
  double * nX = new double[system->nparams()];
  double * dJ = new double[system->nparams()*system->nparams()];

  system->RJ(x, R0, dJ);
  
  for (int i=0; i<system->nparams(); i++) {
    std::copy(x, x+system->nparams(), nX);
    double dx = eps * fabs(nX[i]);
    if (dx < eps) dx = eps;
    nX[i] += dx;
    system->RJ(nX, nR, dJ);
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
double diff_jac_check(std::shared_ptr<Solvable> system, const double * const x,
                      const double * const J)
{
  double * nJ = new double[system->nparams() * system->nparams()];
  
  diff_jac(system, x, nJ);
  double ss = 0.0;
  double js = 0.0;
  for (int i=0; i< system->nparams() * system->nparams(); i++) {
    ss += pow(J[i] - nJ[i], 2.0);
    js += pow(J[i], 2.0);
  }

  delete [] nJ;

  return ss/js;
}


} // namespace neml
