#ifndef NEMLMATH_H
#define NEMLMATH_H

#include <cstddef>

#define CINDEX(i,j,n) (j + i * n)

// Note all this uses the notation:
// symmetric: [E(0,0), E(1,1), E(2,2), sqrt(2) E(1,2), sqrt(2) E(0,2), sqrt(2) E(0,1)]
// skew: [-W(1,2), W(0,2), -W(0,1)]

// BLAS/lapack defs
extern "C" {
  void dgetrf_(const int & m, const int & n, double* A, const int & lda, int* ipiv, int & info);
  void dgetri_(const int & n, double* A, const int & lda, int* ipiv, double* work, const int & lwork, int & info);
  void dgesv_(const int & n, const int & nrhs, double * A, const int & lda, int * ipiv, double * b, const int & ldb, int & info);
  void dgemv_(const char * trans, const int & m, const int & n, const double & alpha, const double * A, const int & lda, const double * x, const int & incx, const double & beta, double * y, const int & incy);
  void dgemm_(const char * transa, const char * transb, const int & m, const int & n, const int & k, const double & alpha, const double * A, const int & lda, const double * B, const int & ldb, const double & beta, double * C, const int & ldc);
  void dger_(const int & m, const int & n, const double & alpha, const double * x, const int & incx, const double * y, const int & incy, double * A, const int & lda);
  void dgecon_(const char * norm, const int & n, double * A, const int & lda, const double & anrom, double & rcond, double * work, int * lwork, int & info);
  void dgttrf_(const int & N, double * DL, double * D, double * DU, double * DU2, int * IPIV, int & INFO);
  void dgttrs_(const char * TRANS, const int & N, const int & NRHS, double * DL, double * D, double * DU, double * DU2, int * IPIV, double * B, const int & LDB, int & info);
}

namespace neml {

/// E_im W_kl - W_ik E_jn
int skew_product(const double * const e, const double * const w, double * const A);

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

/// Return the deviatoric vector
int dev_vec(double * const a);

/// Outer product of two vectors
int outer_vec(const double * const a, int na, const double * const b, int nb, double * const C);

/// Rank 2 update
int outer_update(const double * const a, int na, const double * const b, int nb, double * const C);

/// Rank 2 update
int outer_update_minus(const double * const a, int na, const double * const b, int nb, double * const C);

/// Matrix-vector c = A . b
int mat_vec(const double * const A, int m, const double * const b, int n, double * const c);

/// Matrix-vector c = A.T . b
int mat_vec_trans(const double * const A, int m, const double * const b, int n, double * const c);

// Matrix-matrix C = A . B
int mat_mat(int m, int n, int k, const double * const A, const double * const B, double * const C);

/// Invert a matrix in place
int invert_mat(double* const A, int n);

/// Solve unsymmetric system
int solve_mat(const double * const A, int n, double * const x);

/// Get the condition number of a matrix
double condition(const double * const A, int n);

/// Evaluate a polynomial with Horner's method, highest order term first
double polyval(const double * const poly, const int n, double x);

}

#endif // NEMLMATH_H
