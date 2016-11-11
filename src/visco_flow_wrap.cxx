#include "visco_flow.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {
PYBIND11_PLUGIN(visco_flow) {
  py::module m("visco_flow", "Viscoplastic flow models.");

  py::class_<GFlow, std::shared_ptr<GFlow>>(m, "GFlow")
      .def("g", &GFlow::g, "g function in Perzyna model")
      .def("dg", &GFlow::dg, "Derivative of g wrt f")
      ;

  py::class_<GPowerLaw, std::shared_ptr<GPowerLaw>>(m, "GPowerLaw", py::base<GFlow>())
      .def(py::init<double>(), py::arg("n"))
      .def_property_readonly("n", &GPowerLaw::n)
      ;



  py::class_<ViscoPlasticFlowRule, std::shared_ptr<ViscoPlasticFlowRule>>(m, "RateIndenpendentFlowRule")
      .def_property_readonly("nhist", &ViscoPlasticFlowRule::nhist, "Number of history variables.")
      .def("init_hist",
           [](ViscoPlasticFlowRule & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nhist());
            int ier = m.init_hist(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize history variables.")

      .def("y",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> double
           {
            double fv;
            int ier = m.y(arr2ptr<double>(s), arr2ptr<double>(alpha), T, fv);
            py_error(ier);
            return fv;
           }, "Plastic multiplier value.")
      .def("dy_ds",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            int ier = m.dy_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Plastic multiplier derivative with respect to stress.")
      .def("dy_da",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            int ier = m.dy_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Plastic multiplier derivative with respect to history.")

      .def("g",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            int ier = m.g(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Flow rule.")
      .def("dg_ds",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            int ier = m.dg_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Flow rule derivative with respect to stress.")
      .def("dg_da",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,m.nhist());
            int ier = m.dg_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Flow rule derivative with respect to history.")

      .def("h",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            int ier = m.h(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule.")
      .def("dh_ds",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            int ier = m.dh_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule derivative with respect to stress.")
      .def("dh_da",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            int ier = m.dh_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule derivative with respect to history.")
      ;

    py::class_<PerzynaFlowRule, std::shared_ptr<PerzynaFlowRule>>(m, "PerzynaFlowRule", py::base<ViscoPlasticFlowRule>())
        .def(py::init<std::shared_ptr<YieldSurface>, std::shared_ptr<HardeningRule>, std::shared_ptr<GFlow>, double>(),
             py::arg("surface"), py::arg("hardening"), py::arg("g"), py::arg("eta"))
        .def_property_readonly("eta", &PerzynaFlowRule::eta)
        ;

  return m.ptr();
}

} // namespace neml
