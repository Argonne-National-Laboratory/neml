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
  
  py::class_<Interpolate, std::shared_ptr<Interpolate>>(m, "Interpolate")
      .def("value", &Interpolate::value, "Interpolate to x")
      .def("derivative", &Interpolate::derivative, "Derivative at x")
      .def("__call__", 
           [](Interpolate & m, double x) -> double
           {
            return m(x);
           }, "Operator overload ()")
      .def_property_readonly("valid", &Interpolate::valid)
      ;

  py::class_<InvalidInterpolate, std::shared_ptr<InvalidInterpolate>>(m, "InvalidInterpolate", py::base<Interpolate>())
      .def(py::init<>())
      ;

  py::class_<PolynomialInterpolate, std::shared_ptr<PolynomialInterpolate>>(m, "PolynomialInterpolate", py::base<Interpolate>())
      .def(py::init<std::vector<double>>(), py::arg("coefs"))
      ;

  py::class_<PiecewiseLinearInterpolate, std::shared_ptr<PiecewiseLinearInterpolate>>(m, "PiecewiseLinearInterpolate", py::base<Interpolate>())
      .def(py::init<std::vector<double>, std::vector<double>>(),
           py::arg("points"), py::arg("values"))
      ;

  py::class_<PiecewiseLogLinearInterpolate, std::shared_ptr<PiecewiseLogLinearInterpolate>>(m, "PiecewiseLogLinearInterpolate", py::base<Interpolate>())
      .def(py::init<std::vector<double>, std::vector<double>>(),
           py::arg("points"), py::arg("values"))
      ;

  py::class_<ConstantInterpolate, std::shared_ptr<ConstantInterpolate>>(m, "ConstantInterpolate", py::base<Interpolate>())
      .def(py::init<double>(), py::arg("value"))
      ;

  py::class_<MTSShearInterpolate, std::shared_ptr<MTSShearInterpolate>>(m, "MTSShearInterpolate", py::base<Interpolate>())
      .def(py::init<double, double, double>(), py::arg("y0"), py::arg("D"), 
           py::arg("x0"))
      ;

  return m.ptr();
}



} // namespace neml
