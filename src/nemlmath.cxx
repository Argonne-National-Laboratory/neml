#include "nemlmath.h"

namespace neml {

int invert_matrix(double* const A, int n)
{
  int * ipiv = new int[n + 1];
  int lwork = n * n;
  double * work = new double[lwork];
  int info;

  dgetrf_(&n, &n, A, &n, ipiv, &info);
  dgetri_(&n, A, &n, ipiv, work, &lwork, &info);

  delete [] ipiv;
  delete [] work;

  return 0;
}


} // namespace neml
