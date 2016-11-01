#include "nemlerror.h"

namespace neml {

LinalgError::LinalgError(const char * m) :
    std::runtime_error(m)
{

}


} // namespace neml
