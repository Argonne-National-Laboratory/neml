#ifndef NEMLERROR_H
#define NEMLERROR_H

#include "windows.h"

#include <stdexcept>
#include <string>

namespace neml {

class NEML_EXPORT NEMLError: public std::runtime_error {
 public:
  NEMLError(std::string msg);
  std::string message() const;
  
 private:
  std::string msg_;

};

/// Descriptive exception for problems with BLAS/LAPACK
class NEML_EXPORT LinalgError: public NEMLError {
 public:
  LinalgError(std::string msg);
};

/// Descriptive error for when a nonlinear solve failed
class NEML_EXPORT NonlinearSolverError: public NEMLError {
 public:
  NonlinearSolverError(std::string msg);
};

} // namespace neml

#endif // NEMLERROR
