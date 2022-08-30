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
      .def("s", &GeneralFlowRule::s, "Stress rate")
      .def("ds_ds", &GeneralFlowRule::ds_ds, "Stress rate derivative with respect "
           "to stress.")
      .def("ds_da", &GeneralFlowRule::ds_da, "Stress rate derivative with respect "
           "to history.")
      .def("ds_de", &GeneralFlowRule::ds_de, "Stress rate derivative with respect "
           "to strain rate.")
      .def("a", &GeneralFlowRule::a, "History rate.")
      .def("da_ds", &GeneralFlowRule::da_ds, "History rate derivative with respect "
           "to stress.")
      .def("da_da", &GeneralFlowRule::da_da, "History rate derivative with respect "
           "to history.")
      .def("da_de", &GeneralFlowRule::da_de, "History rate derivative with respect "
           "to strain rate.") 
      .def("work_rate", &GeneralFlowRule::work_rate, "Inelastic work rate.")
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
