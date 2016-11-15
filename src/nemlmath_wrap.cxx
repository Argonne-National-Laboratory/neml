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

  m.def("dev_vec",
        [](py::array_t<double, py::array::c_style> a) -> py::array_t<double>
        {
          if (a.request().ndim != 1) {
            throw LinalgError("The array must be a vector!");
          }
          if (a.request().shape[0] != 6) {
            throw LinalgError("a must be a 6-vector!");
          }
          
          int ier = dev_vec(arr2ptr<double>(a));
          py_error(ier);

          return a;

        }, "Make a vector deviatoric, IN PLACE.");

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

  m.def("outer_update",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b, py::array_t<double, py::array::c_style> C) -> py::array_t<double>
        {
          if ((a.request().ndim != 1)  || (b.request().ndim != 1)) {
            throw LinalgError("The arrays must be vectors!");
          }

          if (C.request().ndim != 2) {
            throw LinalgError("C must be a matrix!");
          }

          if ((C.request().shape[0] != a.request().shape[0]) || (C.request().shape[1] != b.request().shape[0])) {
            throw LinalgError("C is not conformable with axb!");
          }
          
          outer_update(arr2ptr<double>(a), a.request().shape[0], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(C));

          return C;
        }, "Rank 2 update C += a x b.");

  m.def("outer_update_minus",
        [](py::array_t<double, py::array::c_style> a, py::array_t<double, py::array::c_style> b, py::array_t<double, py::array::c_style> C) -> py::array_t<double>
        {
          if ((a.request().ndim != 1)  || (b.request().ndim != 1)) {
            throw LinalgError("The arrays must be vectors!");
          }

          if (C.request().ndim != 2) {
            throw LinalgError("C must be a matrix!");
          }

          if ((C.request().shape[0] != a.request().shape[0]) || (C.request().shape[1] != b.request().shape[0])) {
            throw LinalgError("C is not conformable with axb!");
          }
          
          outer_update_minus(arr2ptr<double>(a), a.request().shape[0], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(C));

          return C;
        }, "Rank 2 update C -= a x b.");

  m.def("mat_vec",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("A must be a matrix!");
          }
          if (b.request().ndim != 1) {
            throw LinalgError("b must be a vector!");
          }
          if (A.request().shape[1] != b.request().shape[0]) {
            throw LinalgError("A and b are not conformable!");
          }

          auto c = alloc_vec<double>(A.request().shape[0]);
          
          mat_vec(arr2ptr<double>(A), A.request().shape[0], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(c));

          return c;

        }, "Matrix-vector product c = A.b.");

  m.def("mat_vec_trans",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("A must be a matrix!");
          }
          if (b.request().ndim != 1) {
            throw LinalgError("b must be a vector!");
          }
          if (A.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("A and b are not conformable!");
          }

          auto c = alloc_vec<double>(A.request().shape[1]);
          
          mat_vec_trans(arr2ptr<double>(A), A.request().shape[1], arr2ptr<double>(b), b.request().shape[0], arr2ptr<double>(c));

          return c;

        }, "Matrix-vector product c = A.T.b.");

  m.def("mat_mat",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> B) -> py::array_t<double>
        {
          if ((A.request().ndim != 2) || (B.request().ndim != 2)) {
            throw LinalgError("A and B must be matrices!");
          }
          if (A.request().shape[1] != B.request().shape[0]) {
            throw LinalgError("A and B must be conformal!");
          }

          auto C = alloc_mat<double>(A.request().shape[0], B.request().shape[1]);
          
          mat_mat(A.request().shape[0], B.request().shape[1], A.request().shape[1],
                  arr2ptr<double>(A), arr2ptr<double>(B), arr2ptr<double>(C));

          return C;
        }, "Matrix-matrix product C = A.B.");

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

   m.def("solve_mat",
        [](py::array_t<double, py::array::c_style> A, py::array_t<double, py::array::c_style> b) -> py::array_t<double>
        {
          if (A.request().ndim != 2) {
            throw LinalgError("A is not a matrix!");
          }
          if (A.request().shape[0] != A.request().shape[1]) {
            throw LinalgError("A is not square!");
          }
          if (b.request().ndim != 1) {
            throw LinalgError("b is not a vector!");
          }
          if (A.request().shape[0] != b.request().shape[0]) {
            throw LinalgError("A and b are not conformable!");
          }

          solve_mat(arr2ptr<double>(A), A.request().shape[0], arr2ptr<double>(b));

          return b;
        }, "Solve Ax=b.");

   m.def("condition",
        [](py::array_t<double, py::array::c_style> A) -> double
        {
          if (A.request().ndim != 2) {
            throw LinalgError("A is not a matrix!");
          }
          if (A.request().shape[0] != A.request().shape[1]) {
            throw LinalgError("A is not square!");
          }

          return condition(arr2ptr<double>(A), A.request().shape[1]);
        }, "Calculate the approximate condition number of A.");
}

} // namespace neml
