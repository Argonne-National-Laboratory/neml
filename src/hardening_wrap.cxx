#include "hardening.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(hardening) {
  py::module m("hardening", "Various hardening rules.");

  py::class_<AssociativeHardening, std::shared_ptr<AssociativeHardening>>(m, "AssociativeHardening")
      .def_property_readonly("nhist", &AssociativeHardening::nhist, "Number of history variables.")

      .def("init_hist",
           [](const AssociativeHardening & m) -> py::array_t<double>
           {
            auto v = alloc_vec<double>(m.nhist());
            m.init_hist(arr2ptr<double>(v));
            return v;
           }, "Initialize history.")
        
      .def("q",
           [](const AssociativeHardening & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto q = alloc_vec<double>(m.nhist());
            m.q(arr2ptr<double>(alpha), T, arr2ptr<double>(q));
            return q;
           }, "Map alpha to q.")

      .def("D",
           [](const AssociativeHardening & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto D = alloc_mat<double>(m.nhist(), m.nhist());
            m.D(arr2ptr<double>(alpha), T, arr2ptr<double>(D));
            return D;
           }, "Generalized plastic modulii.")

      .def("D_inv",
           [](const AssociativeHardening & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto Di = alloc_mat<double>(m.nhist(), m.nhist());
            m.D_inv(arr2ptr<double>(alpha), T, arr2ptr<double>(Di));
            return Di;
           }, "Inverse of generalized plastic modulii.")
      ;
  
  py::class_<IsoJ2LinearAHardening, std::shared_ptr<IsoJ2LinearAHardening>>(m, "IsoJ2LinearAHardening", py::base<AssociativeHardening>())
      .def(py::init<double, double>())
      
      .def_property_readonly("K0", &IsoJ2LinearAHardening::K0)
      .def_property_readonly("Kp", &IsoJ2LinearAHardening::Kp)
      ;

  py::class_<IsoJ2VoceAHardening, std::shared_ptr<IsoJ2VoceAHardening>>(m, "IsoJ2VoceAHardening", py::base<AssociativeHardening>())
      .def(py::init<double, double, double>())
      
      .def_property_readonly("K0", &IsoJ2VoceAHardening::K0)
      .def_property_readonly("Ksat", &IsoJ2VoceAHardening::Ksat)
      .def_property_readonly("delta", &IsoJ2VoceAHardening::delta)
      ;

  return m.ptr();

}


} // namespace neml
