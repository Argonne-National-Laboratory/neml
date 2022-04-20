#ifndef NEMLMATH_H
#define NEMLMATH_H

#include "../windows.h"

#include <cstddef>
#include <vector>
#include <string>
#include <cmath>

#define CINDEX(i,j,n) ((j) + (i) * (n))

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
  void dsyev_(const char * JOBZ, const char * UPLO, const int & N, double * A, const int & LDA, double * W, double * WORK, const int & LWORK, const int & INFO);
}

namespace neml {

const double idsym[81] = {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0};
const double idskew[81] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,-0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,-0.5,0.0,0.0,0.0,-0.5,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,-0.5,0.0,0.0,0.0,-0.5,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.5,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

/// Specialty crystal plasticity operator: M_kmab*W_ml - W_km*M_mlab
NEML_EXPORT void SymSymR4SkewmSkewSymR4SymR4(const double * const M, const double * const W, double * const SS);

/// Specialty crystal plasticity operator: D_km * M_mlab - M_kmab D_ml
NEML_EXPORT void SymSkewR4SymmSkewSymR4SymR4(const double * const D, const double * const M, double * const SS);

/// Specialty operator for the skew part of the tangent: C_ijkb * e_ka - C_ijal * e_bl
NEML_EXPORT void SpecialSymSymR4Sym(const double * const D, const double * const M, double * const SW);

/// Convert the symmetric and skew parts voido a complete fourth order
NEML_EXPORT void transform_fourth(const double * const D, const double * const W, double * const M);

/// The outer product used in constructing the truesdell tangent
NEML_EXPORT void truesdell_tangent_outer(const double * const S, double * const M);

/// Convert a 9x9 to a skew derivative matrix
NEML_EXPORT void full2skew(const double * const A, double * const M);

/// Convert a skew derivative matrix to a 9x9
NEML_EXPORT void skew2full(const double * const M, double * const A);

/// Convert a 9x9 voido a ws matrix
NEML_EXPORT void full2wws(const double * const A, double * const M);

/// Convert a 3x6 ws matrix voido a 9x9
NEML_EXPORT void wws2full(const double * const M, double * const A);

/// Convert a 9x9 to a mandel matrix
NEML_EXPORT void full2mandel(const double * const A, double * const M);

/// Convert a mandel matrix to a full 9x9
NEML_EXPORT void mandel2full(const double * const M, double * const A);

/// Convect a symmetric tensor with a Truesdell rate
NEML_EXPORT void truesdell_update_sym(const double * const D, const double * const W,
                         const double * const Sn, const double * const So,
                         double * const Snp1);

/// Form the 9x9 matrix used in the update
NEML_EXPORT void truesdell_mat(const double * const D, const double * const W,
                  double * const M);

/// Form the RHS of the update, as a symmetric vector
NEML_EXPORT void truesdell_rhs(const double * const D, const double * const W,
                  const double * const Sn, const double * const So,
                  double * const St);

/// Convert a full symmetric rank 2 to a Mandel vector
NEML_EXPORT void sym(const double * const A, double * const v);

/// Convert a symmetric vector to a full matrix
NEML_EXPORT void usym(const double * const v, double * const A);

/// Convert a full skew rank 2 to a skew vector
NEML_EXPORT void skew(const double * const A, double * const v);

/// Convert a skew vector to a full matrix
NEML_EXPORT void uskew(const double * const v, double * const A);

/// Negate a vector in place
NEML_EXPORT void minus_vec(double * const a, int n);

/// Add vectors
NEML_EXPORT void add_vec(const double * const a, const double * const b, int n, double * const c);

/// Subtract vectors
NEML_EXPORT void sub_vec(const double * const a, const double * const b, int n, double * const c);

/// Compute a dot product
NEML_EXPORT double dot_vec(const double * const a, const double * const b, int n);

/// Compute a two norm
NEML_EXPORT double norm2_vec(const double * const a, int n);

/// Normalize a vector in place (2-norm)
NEML_EXPORT void normalize_vec(double * const a, int n);

/// Return the deviatoric vector
NEML_EXPORT void dev_vec(double * const a);

/// Outer product of two vectors
NEML_EXPORT void outer_vec(const double * const a, int na, const double * const b, int nb, double * const C);

/// Rank 2 update
NEML_EXPORT void outer_update(const double * const a, int na, const double * const b, int nb, double * const C);

/// Rank 2 update
NEML_EXPORT void outer_update_minus(const double * const a, int na, const double * const b, int nb, double * const C);

/// Matrix-vector c = A . b
NEML_EXPORT void mat_vec(const double * const A, int m, const double * const b, int n, double * const c);

/// Matrix-vector c = A.T . b
NEML_EXPORT void mat_vec_trans(const double * const A, int m, const double * const b, int n, double * const c);

// Matrix-matrix C = A . B
NEML_EXPORT void mat_mat(int m, int n, int k, const double * const A, const double * const B, double * const C);

// Matrix-matrix C = A . B.T
NEML_EXPORT void mat_mat_ABT(int m, int n, int k, const double * const A, const double * const B, double * const C);

/// Invert a matrix in place
NEML_EXPORT void invert_mat(double* const A, int n);

/// Solve unsymmetric system
NEML_EXPORT void solve_mat(const double * const A, int n, double * const x);

/// Get the condition number of a matrix
NEML_EXPORT double condition(const double * const A, int n);

/// Evaluate a polynomial with Horner's method, highest order term first
NEML_EXPORT double polyval(const std::vector<double> & poly, double x);

/// Construct a polynomial with the given roots
NEML_EXPORT std::vector<double> poly_from_roots(const std::vector<double> & roots);

/// Get the derivative of a polynomial
NEML_EXPORT std::vector<double> differentiate_poly(const std::vector<double> & poly,
                                       int n = 1);

/// The greatest common divisor between two numbers
NEML_EXPORT int gcd(int a, int b);

/// The greatest common divisor between a lot of numbers
NEML_EXPORT int common_gcd(std::vector<int> in);

/// Divide a vector by the collective GCD
NEML_EXPORT std::vector<int> reduce_gcd(std::vector<int> in);

/// Convert to radians
NEML_EXPORT double convert_angle(double a, std::string);

/// Convert from radians
NEML_EXPORT double cast_angle(double a, std::string angles);

/// Vectorized quaternion multiplication
NEML_EXPORT void qmult_vec(const double * const As, const double * const B,
               size_t n, double * const C);

#define RTOL 1.0E-5
#define ATOL 1.0E-8

/// This is only to be used for testing
NEML_EXPORT bool isclose(double a, double b);

/// Perform A * B * A.T
NEML_EXPORT void rotate_matrix(int m, int n, const double * const A,
                  const double * const B, double * C);

/// Factorial
NEML_EXPORT int fact(int n);

/// Factorial as a double + cacheing
NEML_EXPORT double factorial(int n);

/// Get the eigenvalues of a symmetric 3x3 matrix in Mandel notation
NEML_EXPORT void eigenvalues_sym(const double * const s, double * values);

/// Get the eigenvectors of a symmetric 3x3 matrix (row major)
NEML_EXPORT void eigenvectors_sym(const double * const s, double * vectors);

/// First principal invariant
NEML_EXPORT double I1(const double * const s);

/// Second principal invariant
NEML_EXPORT double I2(const double * const s);

}

#endif // NEMLMATH_H
