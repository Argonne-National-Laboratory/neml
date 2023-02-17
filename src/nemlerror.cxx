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

void throw_exception(ExceptionType type)
{
  switch (type) {
    case ExceptionType::NEMLError:
      throw NEMLError("This is an error");
      break;
    case ExceptionType::LinalgError:
      throw LinalgError("This is another error");
      break;
    case ExceptionType::NonlinearSolverError:
      throw NonlinearSolverError("This is still another error");
      break;
  }
}

} // namespace neml
