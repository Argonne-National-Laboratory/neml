#ifndef NEMLMATH_H
#define NEMLMATH_H

#include <cstddef>

#define CINDEX(i,j,n) (j + i * n)

// BLAS/lapack defs
extern "C" {
  void dgetrf_(int* M, int *N, double* A, int* lda, int* ipiv, int* info);
  void dgetri_(int* N, double* A, int* lda, int* ipiv, double* work, int* lwork, int* info);
}

namespace neml {

/// Negate a vector in place
int minus_vec(double * const a, int n);

/// Add vectors
int add_vec(const double * const a, const double * const b, int n, double * const c);

/// Subtract vectors
int sub_vec(const double * const a, const double * const b, int n, double * const c);

/// Compute a dot product
double dot_vec(const double * const a, const double * const b, int n);

/// Compute a two norm
double norm2_vec(const double * const a, int n);

/// Normalize a vector in place (2-norm)
int normalize_vec(double * const a, int n);

/// Outer product of two vectors
int outer_vec(const double * const a, int na, const double * const b, int nb, double * const C);

/// Invert a matrix in place
int invert_mat(double* const A, int n);

}

#endif // NEMLMATH_H
