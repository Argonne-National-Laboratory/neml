#include "harmonics.h"

#include "sphere.h"

const std::complex<double> I{0.0,1.0};

namespace neml {

std::complex<double> harmonic_SO3(int n, int i, int j, const Orientation & pt)
{
  double phi1, Phi, phi2;
  pt.to_euler(phi1, Phi, phi2, "radians", "bunge");

  std::complex<double> a = std::exp(I * (double) i * phi2);
  std::complex<double> b = std::exp(I * (double) j * phi1);
  std::complex<double> P = P_SO3(n, i, j, cos(Phi));

  double nf = sqrt(((double) 2 * n + 1) / (8.0 * M_PI * M_PI));

  return a * P * b * nf;
}

std::complex<double> P_SO3(int n, int i, int j, double x)
{
  // How do people actually handle this?
  if (x == 1.0) {
    if (i == j) {
      return pow(-1.0, (double) n + i);
    }
    else {
      return 0.0;
    }
  }

  std::vector<double> roots;
  for (int ii=0; ii<(n-i); ii++) {
    roots.push_back(1.0);
  }
  for (int ii=0; ii<(n+i); ii++) {
    roots.push_back(-1.0);
  }
  
  std::vector<double> poly = poly_from_roots(roots);

  std::vector<double> pder = differentiate_poly(poly, n - j);

  double val = polyval(pder, x);

  std::complex<double> pf1 = std::pow(-1.0, (double) (n - i)) * 
      std::pow(I, (double) (j - i)) / (std::pow(2.0, (double) n) * factorial(n-i));

  std::complex<double> pf2 = std::sqrt(factorial(n-i)*factorial(n+j) / 
                    (factorial(n+i)*factorial(n-j)));

  std::complex<double> pf3 = std::pow(1.0 - x, (double) -(j-i) / 2.0) * 
      std::pow(1.0 + x, (double) -(j+i) / 2.0);

  return pf1 * pf2 * pf3 * val;
}

void quadrature_S1(size_t n, std::vector<double> & pts, std::vector<double> & wts)
{
  size_t npts = n + 1;
  double di = 2.0 * M_PI / ((double) npts);
  pts.resize(npts);
  wts.resize(npts);
  for (size_t i = 0; i < npts; i++) {
    pts[i] = ((double) i) * di + di / 2;
    wts[i] = di;
  }
}

void quadrature_S2(size_t n, std::vector<double> & theta, 
                   std::vector<double> & phi, std::vector<double> & wts)
{
  size_t np = n;
  if (np > MAX_S2_DEGREE) {
    throw std::invalid_argument("Integration order for S2 too high");
  }
  if (np == 0) {
    theta.resize(1);
    theta[0] = 0.0;
    phi.resize(1);
    phi[0] = 0.0;
  }
  else {
    theta.assign(S2_theta[np-1].begin(), S2_theta[np-1].end());
    phi.assign(S2_phi[np-1].begin(), S2_phi[np-1].end());
  }
  size_t npts = theta.size();
  wts.resize(npts);
  for (size_t i = 0; i<npts; i++) {
    wts[i] = 4.0 * M_PI / ((double) npts);
  }
}

void quadrature_SO3(size_t n, std::vector<Orientation> & pts, 
                    std::vector<double> & wts)
{
  size_t use = n; // Round off causes problems for n>8
  if (use > MAX_SO3_DEGREE) {
    throw std::invalid_argument("Integration order for SO(3) too high");
  }

  for (auto & v : SO3_pts[use]) {
    pts.push_back(Orientation(v));
  }
  
  size_t npts = pts.size();
  wts.resize(npts);
  std::fill(wts.begin(), wts.end(), (8.0 * M_PI * M_PI) / ((double) npts));
}

} // namespace neml
