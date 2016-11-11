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
          double atol, double rtol, int miter, bool verbose)
{
  // For now just plain NR
  return newton(system, x, atol, rtol, miter, verbose);
}

int newton(std::shared_ptr<Solvable> system, double * x,
          double atol, double rtol, int miter, bool verbose)
{
  int n = system->nparams();
  system->init_x(x);
  double * R = new double[n];
  double * J = new double[n*n];
  
  system->RJ(x, R, J);

  double nR = norm2_vec(R, n);
  double nR0 = nR;
  int i = 0;
  int ier = 0;

  if (verbose) {
    std::cout << "Iter.\tnR\t\tnR/nR0" << std::endl;
    std::cout << std::setw(6) << std::left << i << "\t" << std::setw(8) << std::left << std::scientific << nR << "\t" << std::setw(8) << std::left << nR/nR0 << std::endl;
  }

  while ((nR > atol) && (i < miter))
  {
    solve_mat(J, n, R);
    for (int j=0; j<n; j++) x[j] -= R[j];

    system->RJ(x, R, J);
    nR = norm2_vec(R, n);
    i++;

    if (verbose) {
      std::cout << i << "\t" << nR << "\t" << nR/nR0 << std::endl;
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

} // namespace neml
