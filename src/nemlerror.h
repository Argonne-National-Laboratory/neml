#ifndef NEMLERROR_H
#define NEMLERROR_H

#include <stdexcept>

namespace neml {

class LinalgError: public std::runtime_error {
 public:
  LinalgError(const char* m);

};


} // namespace neml

#endif // NEMLERROR
