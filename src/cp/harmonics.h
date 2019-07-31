#ifndef HARMONICS_H
#define HARMONICS_H

#include <complex>

#include "../math/rotations.h"
#include "../math/nemlmath.h"

namespace neml {

/// Generalized spherical harmonic n <= i <= n, -n <= j <= n
std::complex<double> harmonic_SO3(int n, int i, int j, const Orientation & pt);

/// Associated Legendre polynomial for the SO(3) harmonics
std::complex<double> P_SO3(int n, int i, int j, double x);

/// Quadrature rule for S1
void quadrature_S1(size_t n, std::vector<double> & pts, std::vector<double> & wts);

/// Quadrature rule for S2
void quadrature_S2(size_t n, std::vector<double> & theta, 
                   std::vector<double> & phi, std::vector<double> & wts);

/// Quadrature rule for SO(3)
void quadrature_SO3(size_t n, std::vector<Orientation> & pts, 
                    std::vector<double> & wts);

} // namespace neml

#endif // HARMONICS_H
