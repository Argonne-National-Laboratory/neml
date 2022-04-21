#include "nemlerror.h"

#include <iostream>

namespace neml {

NEMLError::NEMLError(std::string msg) :
    std::runtime_error(msg.c_str()), msg_(msg)
{

}

std::string NEMLError::message() const
{
  return msg_;
}

LinalgError::LinalgError(std::string msg) :
    NEMLError(msg)
{

}

NonlinearSolverError::NonlinearSolverError(std::string msg) :
    NEMLError(msg)
{

}


} // namespace neml
