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

/// Invert a matrix in place
int invert_matrix(double* const A, int n);

}

#endif // NEMLMATH_H
