#include "../pyhelp.h"

#include "tensors.h"

#include <iostream>
#include <string>
#include <stdexcept>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(tensors, m) {
  m.doc() = "Standard tensor properties";

  py::class_<Tensor, std::shared_ptr<Tensor>>(m, "Tensor")
      .def(py::self *= double())
      .def(py::self /= double())

      .def(py::self == py::self)
      .def(py::self != py::self)
      ;

  py::class_<Vector, Tensor, std::shared_ptr<Vector>>(m, "Vector")
      .def(py::init<const std::vector<double>>(), py::arg("data"))

      .def_property_readonly("data",
           [](Vector & me) -> py::array_t<double>
           {
            auto v = alloc_vec<double>(3);
            std::copy(me.data(), me.data()+3, arr2ptr<double>(v));
            return v;
           }, "raw data view")

      .def("__repr__",
           [](Vector & me) -> std::string
           {
              std::ostringstream ss;
              ss << "Vector(array([" << me.data()[0] << ", "
                << me.data()[1] << ", " << me.data()[2] << ", " << "]))";
              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](Vector & me) -> std::string
           {
              std::ostringstream ss;
              ss << "[" << me.data()[0] << ", "
                << me.data()[1] << ", " << me.data()[2] << ", " << "]";
              return ss.str();
           }, "python __str__")

      .def("dot", &Vector::dot)
      .def("norm", &Vector::norm)

      .def("opposite", &Vector::opposite)
      .def("__neg__", &Vector::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      .def("cross", &Vector::cross)
      .def("normalize", &Vector::normalize)
      ;

} // PYBIND11_MODULE(tensors, m)

} // namespace neml
