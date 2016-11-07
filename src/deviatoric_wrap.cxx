#include "deviatoric.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include <iostream>

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(deviatoric) {
  py::module m("deviatoric", "Models for the deviatoric stress.");
  
  py::class_<DeviatoricModel, std::shared_ptr<DeviatoricModel>>(m, "DeviatoricModel")
      .def_property_readonly("nhist", &DeviatoricModel::nhist, "Number of history variables.")

      .def("init_hist",
           [](const DeviatoricModel & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nhist());

            int ier = m.init_hist(arr2ptr<double>(h));
            py_error(ier);

            return h;
           }, "Initialize history.")

      .def("update",
           [](const DeviatoricModel & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nhist());
            auto A_np1 = alloc_mat<double>(6,6);
            
            int ier = m.update(arr2ptr<double>(e_np1), arr2ptr<double>(e_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1));
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1);
           }, "Update to next deviatoric stress state.")
      ;

  py::class_<LEModel, std::shared_ptr<LEModel>>(m, "LEModel", py::base<DeviatoricModel>())
      .def(py::init<std::shared_ptr<ShearModulus>>())
      ;

  py::class_<RIAFModel, std::shared_ptr<RIAFModel>>(m, "RIAFModel", py::base<DeviatoricModel>())
      .def(py::init<std::shared_ptr<ShearModulus>, 
            std::shared_ptr<YieldSurface>, 
            std::shared_ptr<AssociativeHardening>,
            double, int, bool>(), py::arg("shear_modulus"),
            py::arg("yield_surface"), py::arg("hardening_rule"),
            py::arg("tol") = 1.0e-6,
            py::arg("miter") = 25, py::arg("verbose") = false)
      ;

  return m.ptr();

}


} // namespace neml
