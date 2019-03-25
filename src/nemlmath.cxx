#include "nemlmath.h"

#include "nemlerror.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace neml {

int skew_product(const double * const e, const double * const w, double * const A)
{
  std::fill(A, A+36, 0.0);
}

int minus_vec(double * const a, int n)
{
  for (int i=0; i<n; i++) {
    a[i] = -a[i];
  }

  return 0;
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
  if (fabs(nv) < std::numeric_limits<double>::epsilon()) {
    std::fill(a, a+n, 0.0);
    return 0;
  }
  for (int i=0; i<n; i++) {
    a[i] /= nv;
  }
  return 0;
}

int dev_vec(double * const a)
{
  double tr = (a[0] + a[1] + a[2]) / 3.0;
  for (int i=0; i<3; i++) {
    a[i] -= tr;
  }

  return 0;
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

// Problem?
int outer_update(const double * const a, int na, const double * const b, 
                 int nb, double * const C)
{
  dger_(nb, na, 1.0, b, 1, a, 1, C, nb);

  return 0;
}

int outer_update_minus(const double * const a, int na, const double * const b, 
                       int nb, double * const C)
{
  dger_(nb, na, -1.0, b, 1, a, 1, C, nb);

  return 0;
}

int mat_vec(const double * const A, int m, const double * const b, int n, 
            double * const c)
{
  dgemv_("T", n, m, 1.0, A, n, b, 1, 0.0, c, 1);

  return 0;
}

int mat_vec_trans(const double * const A, int m, const double * const b, int n, 
            double * const c)
{
  dgemv_("N", m, n, 1.0, A, m, b, 1, 0.0, c, 1);

  return 0;
}

int invert_mat(double * const A, int n)
{
  int * ipiv = new int[n + 1];
  int lwork = n * n;
  double * work = new double[lwork];
  int info;

  dgetrf_(n, n, A, n, ipiv, info);
  if (info > 0) {
    delete [] ipiv;
    delete [] work;
    return LINALG_FAILURE;
  }

  dgetri_(n, A, n, ipiv, work, lwork, info);

  delete [] ipiv;
  delete [] work;

  if (info > 0) return LINALG_FAILURE;

  return 0;
}

int mat_mat(int m, int n, int k, const double * const A,
            const double * const B, double * const C)
{
  dgemm_("N", "N", n, m, k, 1.0, B, n, A, k, 0.0, C, n);

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

  if (info > 0) return LINALG_FAILURE;
  
  return 0;
}

/*
 *  No error checking in this function, as it is assumed to be non-critical
 */
double condition(const double * const A, int n)
{
  // Setup
  int info;
  int * ipiv = new int [n];
  double * x = new double [n];
  double * B = new double [n*n];
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      B[CINDEX(i,j,n)] = A[CINDEX(j,i,n)];
    }
  }
  
  // 1 norm
  double anorm = 0.0;
  double csum;
  for (int i=0; i<n; i++) {
    csum = 0.0;
    for (int j=0; j<n; j++) {
      csum += fabs(A[CINDEX(i,j,n)]);
    }
    if (csum > anorm) anorm = csum;
  }

  // Solve
  std::fill(x, x+n, 0.0);
  dgesv_(n, 1, B, n, ipiv, x, n, info);
  delete [] ipiv;

  double * work = new double [4*n];
  int * iwork = new int [n];
  double rcond;
  dgecon_("1", n, B, n, anorm, rcond, work, iwork, info);

  delete [] B;
  delete [] work;
  delete [] iwork;

  return 1.0 / rcond;
}

double polyval(const double * const poly, const int n, double x)
{
  double res = poly[0];
  for (int i=1; i < n; i++) {
    res = res * x + poly[i];
  }
  return res;
}

} // namespace neml
