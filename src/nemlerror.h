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
  KT_VIOLATION = -4
} Error;

// Translate to exception
void py_error(int ier);

// Translate to string
std::string string_error(int ier);

// Custom exceptions
class LinalgError: public std::runtime_error {
 public:
  LinalgError(const char* m);

};

} // namespace neml

#endif // NEMLERROR
