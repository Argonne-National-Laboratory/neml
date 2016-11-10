#include "ri_flow.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(ri_flow) {
  py::module m("ri_flow", "Rate independent flow models.");
  
  py::class_<RateIndependentFlowRule, std::shared_ptr<RateIndependentFlowRule>>(m, "RateIndenpendentFlowRule")
      .def_property_readonly("nhist", &RateIndependentFlowRule::nhist, "Number of history variables.")
      .def("init_hist",
           [](RateIndependentFlowRule & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nhist());
            int ier = m.init_hist(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize history variables.")

      .def("f",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> double
           {
            double fv;
            int ier = m.f(arr2ptr<double>(s), arr2ptr<double>(alpha), T, fv);
            py_error(ier);
            return fv;
           }, "Yield surface value.")
      .def("df_ds",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            int ier = m.df_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Yield surface derivative with respect to stress.")
      .def("df_da",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            int ier = m.df_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Yield surface derivative with respect to history.")

      .def("g",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            int ier = m.g(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Flow rule.")
      .def("dg_ds",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            int ier = m.dg_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Flow rule derivative with respect to stress.")
      .def("dg_da",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,m.nhist());
            int ier = m.dg_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Flow rule derivative with respect to history.")

      .def("h",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            int ier = m.h(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule.")
      .def("dh_ds",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            int ier = m.dh_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule derivative with respect to stress.")
      .def("dh_da",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            int ier = m.dh_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule derivative with respect to history.")
      ;

  py::class_<RateIndependentAssociativeFlow, std::shared_ptr<RateIndependentAssociativeFlow>>(m, "RateIndependentAssociativeFlow", py::base<RateIndependentFlowRule>())
      .def(py::init<std::shared_ptr<YieldSurface>, std::shared_ptr<HardeningRule>>(), 
           py::arg("surface"), py::arg("hardening"))
      ;

  py::class_<RateIndependentNonAssociativeHardening, std::shared_ptr<RateIndependentNonAssociativeHardening>>(m, "RateIndependentNonAssociativeHardening", py::base<RateIndependentFlowRule>())
      .def(py::init<std::shared_ptr<YieldSurface>, std::shared_ptr<NonAssociativeHardening>>(), 
           py::arg("surface"), py::arg("hardening"))
      ;

  return m.ptr();
}

} // namespace neml
