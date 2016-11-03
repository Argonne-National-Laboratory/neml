#include "nemlmath.h"

#include "nemlerror.h"
#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(nemlmath) {
  py::module m("nemlmath", "Various helper math functions.");

  m.def("minus_vec",
        [](py::array_t<double, py::array::c_style> a) -> py::array_t<double>
        {
          if (a.request().ndim != 1) {
            throw LinalgError("Array must be a vector!");
          }
          minus_vec(arr2ptr<double>(a), a.request().shape[0]);

          return a;
        }, "Negate a vector IN PLACE.");

  m.def("add_vec",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if ((a.request().ndim != 1) || (b.request().ndim != 1)) {
            throw LinalgError("Both arrays must be vectors!");
          }

          if (a.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("Both vectors must have the same length!");
          }

          auto c = alloc_vec<double>(a.request().shape[0]);

          add_vec(arr2ptr<double>(a), arr2ptr<double>(b), a.request().shape[0], arr2ptr<double>(c));

          return c;

        }, "Add two vectors.");

  m.def("sub_vec",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if ((a.request().ndim != 1) || (b.request().ndim != 1)) {
            throw LinalgError("Both arrays must be vectors!");
          }

          if (a.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("Both vectors must have the same length!");
          }

          auto c = alloc_vec<double>(a.request().shape[0]);

          sub_vec(arr2ptr<double>(a), arr2ptr<double>(b), a.request().shape[0], arr2ptr<double>(c));

          return c;

        }, "Subtract two vectors.");

  m.def("dot_vec",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b) -> double
        {
          if ((a.request().ndim != 1) || (b.request().ndim != 1)) {
            throw LinalgError("Both arrays must be vectors!");
          }

          if (a.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("Both vectors must have the same length!");
          }

          return dot_vec(arr2ptr<double>(a), arr2ptr<double>(b), a.request().shape[0]);

        }, "Compute a dot product between two vectors.");

  m.def("norm2_vec",
        [](py::array_t<double, py::array::c_style> a) -> double
        {
          if (a.request().ndim != 1) {
            throw LinalgError("The array must be a vector!");
          }

          return norm2_vec(arr2ptr<double>(a), a.request().shape[0]);
        }, "Compute the two norm of a vector.");

  m.def("normalize_vec",
        [](py::array_t<double, py::array::c_style> a) -> py::array_t<double>
        {
          if (a.request().ndim != 1) {
            throw LinalgError("The array must be a vector!");
          }
          normalize_vec(arr2ptr<double>(a), a.request().shape[0]);

          return a;
        }, "Normalize a vector IN PLACE.");

  m.def("outer_vec",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if ((a.request().ndim != 1)  || (b.request().ndim != 1)) {
            throw LinalgError("The arrays must be vectors!");
          }
          
          auto C = alloc_mat<double>(a.request().shape[0], b.request().shape[0]);

          outer_vec(arr2ptr<double>(a), a.request().shape[0], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(C));

          return C;
        }, "Outer product of a x b.");

  m.def("invert_mat", 
        [](py::array_t<double, py::array::c_style> A) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("Array is not a matrix!");
          }
          if (A.request().shape[0] != A.request().shape[1]) {
            throw LinalgError("Matrix is not square!");
          }
          int ier = invert_mat(arr2ptr<double>(A), A.request().shape[0]);
          // Should check non-singular
          return A;
        }, "Invert a matrix IN PLACE.");
}

} // namespace neml
