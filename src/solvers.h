#ifndef SOLVERS_H
#define SOLVERS_H

#include <cstddef>
#include <memory>

namespace neml {

/// Generic nonlinear solver interface
class Solvable {
 public:
  virtual size_t nparams() const = 0;
  virtual int init_x(double * const x) = 0;
  virtual int RJ(const double * const x, double * const R, double * const J) = 0;
};

// This is entirely for testing
class TestRosenbrock: public Solvable {
 public:
  TestRosenbrock(size_t N);

  virtual size_t nparams() const;
  virtual int init_x(double * const x);
  virtual int RJ(const double * const x, double * const R, double * const J);

 private:
  const size_t N_;
};

/// Call the built-in solver
int solve(std::shared_ptr<Solvable> system, double * x, 
          double atol = 1.0e-10, double rtol = 1.0e-6, int miter = 25,
          bool verbose = false);

/// Default solver: plain NR
int newton(std::shared_ptr<Solvable> system, double * x, 
          double atol, double rtol, int miter, bool verbose);

} // namespace neml

#endif // SOLVERS_H
