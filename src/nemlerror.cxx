#include "nemlerror.h"

#include <iostream>

namespace neml {

void py_error(int ier)
{
  if (ier != 0)
    throw NEMLError("Unknown error!");
}

std::string string_error(int ier)
{
  return "Unknown error";
}

NEMLError::NEMLError(std::string msg) :
    std::runtime_error(msg.c_str())
{

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
