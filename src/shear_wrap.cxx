#include "shear.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include <memory>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(shear) {
  py::module m("shear", "Shear modulus models, as functions of temperature.");

  py::class_<ShearModulus, std::shared_ptr<ShearModulus>>(m, "ShearModulus")
      .def("modulus", &ShearModulus::modulus, "Return the current shear modulus as a function of temperature.")
  ;

  py::class_<ConstantShearModulus, std::shared_ptr<ConstantShearModulus>>(m, "ConstantShearModulus", py::base<ShearModulus>())
      .def(py::init<double>())
      .def_property_readonly("mu", &ConstantShearModulus::mu)
  ;

}

} // namespace neml
