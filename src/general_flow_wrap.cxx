#include "pyhelp.h" // include first to avoid annoying redef warning

#include "general_flow.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {
PYBIND11_MODULE(general_flow, m) {
  py::module::import("neml.objects");

  m.doc() = "General flow models where subclass functions define everything.";

  py::class_<GeneralFlowRule, HistoryNEMLObject, std::shared_ptr<GeneralFlowRule>>(m, "GeneralFlowRule")
      .def("s",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            m.s(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot, 
                          arr2ptr<double>(f));
            return f;
           }, "Stress rate.")
      .def("ds_ds",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            m.ds_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot, arr2ptr<double>(f));
            return f;
           }, "Stress rate derivative with respect to stress.")
      .def("ds_da",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,m.nhist());
            m.ds_da(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  arr2ptr<double>(f));
            return f;
           }, "Stress rate derivative with respect to history.")
      .def("ds_de",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            m.ds_de(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  arr2ptr<double>(f));
            return f;
           }, "Stress rate derivative with respect to strain.")


      .def("a",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.a(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot, 
                          arr2ptr<double>(f));
            return f;
           }, "History rate.")
      .def("da_ds",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.da_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot, arr2ptr<double>(f));
            return f;
           }, "History rate derivative with respect to stress.")
      .def("da_da",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            m.da_da(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  arr2ptr<double>(f));
            return f;
           }, "History rate derivative with respect to history.")
      .def("da_de",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.da_de(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  arr2ptr<double>(f));
            return f;
           }, "History rate derivative with respect to strain.")
      
      .def("work_rate",
           [](GeneralFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, py::array_t<double, py::array::c_style> edot, double T, double Tdot) -> double
           {
            double pi;
            m.work_rate(arr2ptr<double>(s), arr2ptr<double>(alpha), 
                          arr2ptr<double>(edot), T,
                          Tdot,  pi);
            return pi;
           }, "Plastic work rate.")
      .def("set_elastic_model", &GeneralFlowRule::set_elastic_model)
  ;

  py::class_<TVPFlowRule, GeneralFlowRule, std::shared_ptr<TVPFlowRule>>(m, "TVPFlowRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<TVPFlowRule>(args, kwargs, {"elastic", "flow"});
        }))
      ;
}

} // namespace neml
