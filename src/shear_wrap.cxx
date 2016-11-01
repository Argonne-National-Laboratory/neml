#include "shear.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

namespace neml {

PYBIND11_PLUGIN(shear) {
  py::module m("shear", "Shear modulus models, as functions of temperature.");

  py::class_<ShearModulus>(m, "ShearModulus")
      .def("modulus", &ShearModulus::modulus, "Return the current shear modulus as a function of temperature.")
  ;

  py::class_<ConstantShearModulus>(m, "ConstantShearModulus", py::base<ShearModulus>())
      .def(py::init<double>())
      .def_property_readonly("mu", &ConstantShearModulus::mu)
  ;

}

} // namespace neml
