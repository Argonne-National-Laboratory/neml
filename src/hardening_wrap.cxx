#include "hardening.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(hardening) {
  py::module m("hardening", "Various hardening rules.");

  py::class_<HardeningRule, std::shared_ptr<HardeningRule>>(m, "HardeningRule")
      .def_property_readonly("nhist", &HardeningRule::nhist, "Number of history variables.")

      .def("init_hist",
           [](const HardeningRule & m) -> py::array_t<double>
           {
            auto v = alloc_vec<double>(m.nhist());
            m.init_hist(arr2ptr<double>(v));
            return v;
           }, "Initialize history.")
        
      .def("q",
           [](const HardeningRule & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto q = alloc_vec<double>(m.nhist());
            m.q(arr2ptr<double>(alpha), T, arr2ptr<double>(q));
            return q;
           }, "Map alpha to q.")

      .def("dq_da",
           [](const HardeningRule & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto D = alloc_mat<double>(m.nhist(), m.nhist());
            m.dq_da(arr2ptr<double>(alpha), T, arr2ptr<double>(D));
            return D;
           }, "Gradient of map")
      ;
  
  py::class_<IsotropicHardeningRule, std::shared_ptr<IsotropicHardeningRule>>(m, "IsotropicHardeningRule", py::base<HardeningRule>())
      ;

  py::class_<LinearIsotropicHardeningRule, std::shared_ptr<LinearIsotropicHardeningRule>>(m, "LinearIsotropicHardeningRule", py::base<IsotropicHardeningRule>())
      .def(py::init<double, double>(), py::arg("s0"), py::arg("K"))
      
      .def_property_readonly("s0", &LinearIsotropicHardeningRule::s0)
      .def_property_readonly("K", &LinearIsotropicHardeningRule::K)
      ;

  py::class_<VoceIsotropicHardeningRule, std::shared_ptr<VoceIsotropicHardeningRule>>(m, "VoceIsotropicHardeningRule", py::base<IsotropicHardeningRule>())
      .def(py::init<double, double, double>(), py::arg("s0"), py::arg("R"), py::arg("d"))
      
      .def_property_readonly("s0", &VoceIsotropicHardeningRule::s0)
      .def_property_readonly("R", &VoceIsotropicHardeningRule::R)
      .def_property_readonly("d", &VoceIsotropicHardeningRule::d)
      ;

  py::class_<KinematicHardeningRule, std::shared_ptr<KinematicHardeningRule>>(m, "KinematicHardeningRule", py::base<HardeningRule>())
      ;

  py::class_<LinearKinematicHardeningRule, std::shared_ptr<LinearKinematicHardeningRule>>(m, "LinearKinematicHardeningRule", py::base<KinematicHardeningRule>())
      .def(py::init<double>(), py::arg("H"))
      
      .def_property_readonly("H", &LinearKinematicHardeningRule::H)
      ;

  py::class_<CombinedHardeningRule, std::shared_ptr<CombinedHardeningRule>>(m, "CombinedHardeningRule", py::base<HardeningRule>())
      .def(py::init<std::shared_ptr<IsotropicHardeningRule>, std::shared_ptr<KinematicHardeningRule>>(),
           py::arg("iso"), py::arg("kin"))
      ;

  return m.ptr();

}


} // namespace neml
