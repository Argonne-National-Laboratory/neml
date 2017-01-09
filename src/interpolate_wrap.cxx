#include "interpolate.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(interpolate) {
  py::module m("interpolate", "Interpolation schemes for model parameters.");
  
  py::class_<Interpolate>(m, "Interpolate")
      .def("value", &Interpolate::value, "Interpolate to x")
      .def("__call__", 
           [](Interpolate & m, double x) -> double
           {
            return m(x);
           }, "Operator overload ()")
      ;

  py::class_<PolynomialInterpolate>(m, "PolynomialInterpolate", py::base<Interpolate>())
      .def(py::init<std::vector<double>>(), py::arg("coefs"))
      ;

  py::class_<ConstantInterpolate>(m, "ConstantInterpolate", py::base<Interpolate>())
      .def(py::init<double>(), py::arg("value"))
      ;

  return m.ptr();
}



} // namespace neml
