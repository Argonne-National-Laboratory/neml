#include "pyhelp.h" // include first to avoid annoying redef warning

#include "ri_flow.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(ri_flow, m) {
  py::module::import("neml.objects");

  m.doc() = "Rate independent flow models.";
  
  py::class_<RateIndependentFlowRule, HistoryNEMLObject, std::shared_ptr<RateIndependentFlowRule>>(m, "RateIndenpendentFlowRule")
      .def("f",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> double
           {
            double fv;
            m.f(arr2ptr<double>(s), arr2ptr<double>(alpha), T, fv);
            return fv;
           }, "Yield surface value.")
      .def("df_ds",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            m.df_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Yield surface derivative with respect to stress.")
      .def("df_da",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.df_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Yield surface derivative with respect to history.")

      .def("g",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            m.g(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule.")
      .def("dg_ds",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            m.dg_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule derivative with respect to stress.")
      .def("dg_da",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,m.nhist());
            m.dg_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule derivative with respect to history.")

      .def("h",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.h(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule.")
      .def("dh_ds",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.dh_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule derivative with respect to stress.")
      .def("dh_da",
           [](RateIndependentFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            m.dh_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule derivative with respect to history.")
      ;

  py::class_<RateIndependentAssociativeFlow, RateIndependentFlowRule, std::shared_ptr<RateIndependentAssociativeFlow>>(m, "RateIndependentAssociativeFlow")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<RateIndependentAssociativeFlow>(args, kwargs, {"surface", "hardening"});
        }))
      ;

  py::class_<RateIndependentNonAssociativeHardening, RateIndependentFlowRule, std::shared_ptr<RateIndependentNonAssociativeHardening>>(m, "RateIndependentNonAssociativeHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<RateIndependentNonAssociativeHardening>(args, kwargs, {"surface", "hardening"});
        }))
      ;
}

} // namespace neml
