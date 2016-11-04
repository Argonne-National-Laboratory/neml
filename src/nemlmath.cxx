#include "nemlmath.h"

#include <cmath>
#include <iostream>

namespace neml {

int minus_vec(double * const a, int n)
{
  for (int i=0; i<n; i++) {
    a[i] = -a[i];
  }
}

int add_vec(const double * const a, const double * const b, int n, double * const c)
{
  for (int i=0; i<n; i++) {
    c[i] = a[i] + b[i];
  }

  return 0;
}

int sub_vec(const double * const a, const double * const b, int n, double * const c)
{
  for (int i=0; i<n; i++) {
    c[i] = a[i] - b[i];
  }
  
  return 0;
}

double dot_vec(const double * const a, const double * const b, int n)
{
  double sum = 0.0;
  for (int i=0; i<n; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

double norm2_vec(const double * const a, int n)
{
  return sqrt(dot_vec(a, a, n));
}

int normalize_vec(double * const a, int n)
{
  double nv = norm2_vec(a, n);
  for (int i=0; i<n; i++) {
    a[i] /= nv;
  }
}

int outer_vec(const double * const a, int na, const double * const b, int nb, double * const C)
{
  for (int i=0; i < na; i++) {
    for (int j=0; j < nb; j++) {
      C[CINDEX(i,j,nb)] = a[i] * b[j];
    }
  }

  return 0;
}

int invert_mat(double * const A, int n)
{
  int * ipiv = new int[n + 1];
  int lwork = n * n;
  double * work = new double[lwork];
  int info;

  dgetrf_(n, n, A, n, ipiv, info);
  dgetri_(n, A, n, ipiv, work, lwork, info);

  delete [] ipiv;
  delete [] work;

  return 0;
}

int factor_sym_mat(double * const A, int n)
{
  int info;

  dpotrf_("L", n, A, n, info);

  return 0;
}

int backsolve_sym_mat(double * const A, int n, double * const x)
{
  int info;

  dpotrs_("L", n, 1, A, n, x, n, info);

  return 0;
}

int solve_mat(const double * const A, int n, double * const x)
{
  int info;
  int * ipiv = new int [n];
  double * B = new double [n*n];
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      B[CINDEX(i,j,n)] = A[CINDEX(j,i,n)];
    }
  }
  
  dgesv_(n, 1, B, n, ipiv, x, n, info);

  delete [] ipiv;
  delete [] B;

}


} // namespace neml
