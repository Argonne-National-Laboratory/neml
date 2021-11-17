#include "pyhelp.h"

#include "math/matrix.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(matrix, m) {
  py::module::import("neml.objects");

  m.doc() = "Generic matrix and vector classes";

  py::class_<FlatVector, std::shared_ptr<FlatVector>>(m, "FlatVector",
                                                      py::buffer_protocol())
      .def(py::init<const std::vector<double>>(), py::arg("data"))
      .def_property_readonly("n", &FlatVector::n)
      .def_property_readonly("owns_data", &FlatVector::owns_data)
      .def_buffer(
          [](FlatVector & m) -> py::buffer_info
          {
            return py::buffer_info(
                m.data(),
                sizeof(double),
                py::format_descriptor<double>::format(),
                1,
                {m.n()},
                {sizeof(double)});
          })
      ;

  py::class_<Matrix, std::shared_ptr<Matrix>>(m, "Matrix",
                                              py::buffer_protocol())
      .def(py::init<size_t, size_t>(), py::arg("m"), py::arg("n"))
      .def_property_readonly("m", &Matrix::m)
      .def_property_readonly("n", &Matrix::n)
      .def_property_readonly("size", &Matrix::size)
      .def("dot", &Matrix::dot)
      .def("__call__", 
           [](Matrix & m, size_t i, size_t j) -> double& {
            return m(i,j);
           })
      .def_buffer(
          [](Matrix & m) -> py::buffer_info
          {
            return py::buffer_info(
                m.data(),
                sizeof(double),
                py::format_descriptor<double>::format(),
                2,
                {m.m(), m.n()},
                {sizeof(double) * m.n(), sizeof(double)});
          })
      ;

  py::class_<SquareMatrix, Matrix, NEMLObject, std::shared_ptr<SquareMatrix>>(m,
                                                                              "SquareMatrix",
                                                                              py::buffer_protocol())
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<SquareMatrix>(args,
                                                                kwargs,
                                                                {"m"});
                    }))
      ;

} // PYBIND11_MODULE(tensors, m)

} // namespace neml
