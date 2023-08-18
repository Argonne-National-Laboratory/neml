#include "pyhelp.h" // include first to avoid annoying redef warning

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(nemlerror, m) {
  m.doc() = "NEML-specific errors and exceptions";
  
  py::register_exception<NEMLError>(m, "NEMLError");
  py::register_exception<LinalgError>(m, "LinalgError");
  py::register_exception<NonlinearSolverError>(m, "NonlinearSolverError");
  
  py::enum_<ExceptionType>(m, "ExceptionType")
      .value("NEMLError", ExceptionType::NEMLError)
      .value("LinalgError", ExceptionType::LinalgError)
      .value("NonlinearSolverError", ExceptionType::NonlinearSolverError);
}

} // namespace neml
