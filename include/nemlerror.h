#ifndef NEMLERROR_H
#define NEMLERROR_H

#include "windows.h"

#include <stdexcept>
#include <string>

namespace neml {

// The global error codes
typedef enum Error {
  SUCCESS = 0,
} Error;

/// Translate an error code to an exception
NEML_EXPORT void py_error(int ier);

/// Translate an error code to a string
NEML_EXPORT std::string string_error(int ier);

class NEML_EXPORT NEMLError: public std::runtime_error {
 public:
  NEMLError(std::string msg);

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
