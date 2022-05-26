#include "pyhelp.h"

#include "cp/inelasticity.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(inelasticity, m) {
  m.doc() = "Inelastic models for crystal plasticity";

  py::class_<InelasticModel, HistoryNEMLObject, std::shared_ptr<InelasticModel>>(m,
                                                              "InelasticModel")
      .def("strength", &InelasticModel::strength)
      .def("d_p", &InelasticModel::d_p)
      .def("d_d_p_d_stress", &InelasticModel::d_d_p_d_stress)
      .def("d_d_p_d_history", &InelasticModel::d_d_p_d_history)
      .def("history_rate", &InelasticModel::history_rate)
      .def("d_history_rate_d_stress", &InelasticModel::d_history_rate_d_stress)
      .def("d_history_rate_d_history", &InelasticModel::d_history_rate_d_history)
      .def("w_p", &InelasticModel::w_p)
      .def("d_w_p_d_stress", &InelasticModel::d_w_p_d_stress)
      .def("d_w_p_d_history", &InelasticModel::d_w_p_d_history)
      .def_property_readonly("use_nye", &InelasticModel::use_nye)
      ;

  py::class_<NoInelasticity, InelasticModel, std::shared_ptr<NoInelasticity>>(m,
                                                                              "NoInelasticity")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<NoInelasticity>(args, kwargs, {});
                    }))
      ;

  py::class_<AsaroInelasticity, InelasticModel, std::shared_ptr<AsaroInelasticity>>(m,
                                                                              "AsaroInelasticity")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<AsaroInelasticity>(args, kwargs, {"rule"});
                    }))
      ;

  py::class_<PowerLawInelasticity, InelasticModel,
      std::shared_ptr<PowerLawInelasticity>>(m, "PowerLawInelasticity")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<PowerLawInelasticity>(args, kwargs, {"A", "n"});
                    }))
      ;

  py::class_<CombinedInelasticity, InelasticModel, std::shared_ptr<CombinedInelasticity>>(m, "CombinedInelasticity")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<CombinedInelasticity>(args, kwargs, {"models"});
                    }))
      ;

} // PYBIND11_MODULE

} // namespace neml
