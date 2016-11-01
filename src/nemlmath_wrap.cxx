#include "nemlmath.h"

#include "nemlerror.h"
#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

namespace neml {

PYBIND11_PLUGIN(nemlmath) {
  py::module m("nemlmath", "Various helper math functions.");

  m.def("invert_matrix", 
        [](py::array_t<double, py::array::c_style> A) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("Array is not a matrix!");
          }
          if (A.request().shape[0] != A.request().shape[1]) {
            throw LinalgError("Matrix is not square!");
          }
          int ier = invert_matrix(arr2ptr<double>(A), A.request().shape[0]);
          // Should check non-singular
          return A;
        }, "Invert a matrix IN PLACE.");
}

} // namespace neml
