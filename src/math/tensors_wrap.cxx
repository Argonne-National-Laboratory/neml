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
      .def_property_readonly("data",
           [](Tensor & me) -> py::array_t<double>
           {
            auto v = alloc_vec<double>(me.n());
            std::copy(me.data(), me.data()+me.n(), arr2ptr<double>(v));
            return v;
           }, "raw data view")

      .def(py::self *= double())
      .def(py::self /= double())

      .def(py::self == py::self)
      .def(py::self != py::self)
      ;

  py::class_<Vector, Tensor, std::shared_ptr<Vector>>(m, "Vector")
      // Start standard
      .def(py::init<const std::vector<double>>(), py::arg("data"))

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


      .def("opposite", &Vector::opposite)
      .def("__neg__", &Vector::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      .def("__getitem__", [](const Vector & v, size_t i) {
           if (i >= v.n()) throw py::index_error();
           return v(i);
      })

      .def("__setitem__", [](Vector & v, size_t i, double val) {
           if (i >= v.n()) throw py::index_error();
           v(i) = val;
      })

      // End standard

      .def("dot", &Vector::dot)
      .def("outer", &Vector::outer)
      .def("norm", &Vector::norm)
      .def("cross", &Vector::cross)
      .def("normalize", &Vector::normalize)
      ;

  m.def("outer", &outer);

  py::class_<RankTwo, Tensor, std::shared_ptr<RankTwo>>(m, "RankTwo")
      // Start standard
      .def(py::init<const std::vector<const std::vector<double>>>(), py::arg("data"))

      .def("__repr__",
           [](RankTwo & me) -> std::string
           {
              std::ostringstream ss;
              
              ss << "RankTwo(array([";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << me.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]))";

              return ss.str();
           }, "python __repr__")

      .def("__str__",
           [](RankTwo & me) -> std::string
           {
              std::ostringstream ss;
              
              ss << "[";
              for (size_t i=0; i<3; i++) {
                ss << "[";

                for (size_t j=0; j<3; j++) {
                  ss << me.data()[i*3+j] << " ";
                }
                ss << "]" << std::endl;
              }
              ss << "]";

              return ss.str();
           }, "python __str__")

      .def("opposite", &RankTwo::opposite)
      .def("__neg__", &RankTwo::opposite)

      .def(double() * py::self)
      .def(py::self * double())

      .def(py::self / double())

      .def(py::self += py::self)
      .def(py::self + py::self)

      .def(py::self -= py::self)
      .def(py::self - py::self)

      // End standard
      ;

} // PYBIND11_MODULE(tensors, m)

} // namespace neml
