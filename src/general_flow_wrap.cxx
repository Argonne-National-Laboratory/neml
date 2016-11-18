#include "general_flow.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {
PYBIND11_PLUGIN(general_flow) {
  py::module m("general_flow", "General flow models.");

  py::class_<GeneralFlowRule, std::shared_ptr<GeneralFlowRule>>(m, "GeneralFlowRule")
      .def_property_readonly("nhist", &GeneralFlowRule::nhist, "Number of history variables.")
      .def("init_hist",
           [](GeneralFlowRule & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nhist());
            int ier = m.init_hist(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize history variables.")

      .def("s",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            int ier = m.s(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot, 
                          arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Stress rate.")
      .def("ds_ds",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            int ier = m.ds_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Stress rate derivative with respect to stress.")
      .def("ds_da",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,m.nhist());
            int ier = m.ds_da(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Stress rate derivative with respect to history.")
      .def("ds_de",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            int ier = m.ds_de(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Stress rate derivative with respect to strain.")


      .def("a",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            int ier = m.a(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot, 
                          arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "History rate.")
      .def("da_ds",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            int ier = m.da_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "History rate derivative with respect to stress.")
      .def("da_da",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            int ier = m.da_da(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "History rate derivative with respect to history.")
      .def("da_de",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            int ier = m.da_de(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "History rate derivative with respect to strain.")
  ;

  py::class_<TVPFlowRule, std::shared_ptr<TVPFlowRule>>(m, "TVPFlowRule", py::base<GeneralFlowRule>())
      .def(py::init<std::shared_ptr<LinearElasticModel>, std::shared_ptr<ViscoPlasticFlowRule>>(),
           py::arg("elastic"), py::arg("flow"))
      ;

  return m.ptr();
}

} // namespace neml
