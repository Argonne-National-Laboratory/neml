#ifndef NEMLERROR_H
#define NEMLERROR_H

#include <stdexcept>
#include <string>

namespace neml {

// The global error codes
typedef enum Error {
  SUCCESS = 0,
  INCOMPATIBLE_MODELS = -1,
  LINALG_FAILURE = -2,
  MAX_ITERATIONS = -3,
  KT_VIOLATION = -4,
  NODE_NOT_FOUND = -5,
  TOO_MANY_NODES = -6,
  ATTRIBUTE_NOT_FOUND = -7,
  UNKNOWN_TYPE = -8,
  BAD_TEXT = -9,
  INVALID_TYPE = -10,
  FILE_NOT_FOUND = -11,
  CREEP_PLASTICITY = -12,
  UNKNOWN_ERROR = -13,
  INCOMPATIBLE_KM = -14,
  DUMMY_ELASTIC = -15,
  INCOMPATIBLE_VECTORS = -16
} Error;

/// Translate an error code to an exception
void py_error(int ier);

/// Translate an error code to a string
std::string string_error(int ier);

/// Descriptive exception for problems with BLAS/LAPACK
class LinalgError: public std::runtime_error {
 public:
  LinalgError(const char* m);

};

} // namespace neml

#endif // NEMLERROR
