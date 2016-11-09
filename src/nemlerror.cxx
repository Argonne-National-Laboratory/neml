#include "nemlerror.h"

#include <iostream>

namespace neml {

void py_error(int ier)
{
  switch(ier) {
    case SUCCESS: return;
    case INCOMPATIBLE_MODELS: throw std::runtime_error("Incompatible submodels");
    case LINALG_FAILURE: throw LinalgError("Generic linear algebra failure");
    case MAX_ITERATIONS: throw std::runtime_error("Maximum iteration count exceeded");
    case KT_VIOLATION: throw std::runtime_error("Integration of rate-independent model resulted in a violation of the Kuhn-Tucker conditions");

    default: throw std::runtime_error("Unknown error!");
  }
}

std::string string_error(int ier)
{
  switch(ier) {
    case SUCCESS: return "Success";
    case INCOMPATIBLE_MODELS: return "Incompatible submodels";
    case LINALG_FAILURE: return "Linear algebra call failed";
    case MAX_ITERATIONS: return "Maximum iteration count exceeded";
    case KT_VIOLATION: return "Integration of rate-independent model resulted in a violation of the Kuhn-Tucker conditions";

    default: return "Unknown error";
  }
}

LinalgError::LinalgError(const char * m) :
    std::runtime_error(m)
{

}


} // namespace neml
