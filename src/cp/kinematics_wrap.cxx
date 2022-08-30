#include "pyhelp.h"

#include "cp/kinematics.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(kinematics, m) {
  m.doc() = "Kinematic models for crystal plasticity";
  
  py::class_<KinematicModel, HistoryNEMLObject, std::shared_ptr<KinematicModel>>(m,
                                                              "KinematicModel")

      .def("strength", &KinematicModel::strength)

      .def("decouple", &KinematicModel::decouple)

      .def("stress_rate", &KinematicModel::stress_rate)
      .def("d_stress_rate_d_stress", &KinematicModel::d_stress_rate_d_stress)
      .def("d_stress_rate_d_d", &KinematicModel::d_stress_rate_d_d)
      .def("d_stress_rate_d_w", &KinematicModel::d_stress_rate_d_w)
      .def("d_stress_rate_d_history", &KinematicModel::d_stress_rate_d_history)

      .def("history_rate", &KinematicModel::history_rate)
      .def("d_history_rate_d_stress", &KinematicModel::d_history_rate_d_stress)
      .def("d_history_rate_d_d", &KinematicModel::d_history_rate_d_d)
      .def("d_history_rate_d_w", &KinematicModel::d_history_rate_d_w)
      .def("d_history_rate_d_history", &KinematicModel::d_history_rate_d_history)

      .def("d_stress_rate_d_d_decouple", &KinematicModel::d_stress_rate_d_d_decouple)
      .def("d_stress_rate_d_w_decouple", &KinematicModel::d_stress_rate_d_w_decouple)

      .def("d_history_rate_d_d_decouple", &KinematicModel::d_history_rate_d_d_decouple)
      .def("d_history_rate_d_w_decouple", &KinematicModel::d_history_rate_d_w_decouple)

      .def("spin", &KinematicModel::spin)

      .def("elastic_strains", &KinematicModel::elastic_strains)

      .def_property_readonly("use_nye", &KinematicModel::use_nye)
      ;

  py::class_<StandardKinematicModel, KinematicModel, std::shared_ptr<StandardKinematicModel>>(m, "StandardKinematicModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
           {
            return create_object_python<StandardKinematicModel>(args, kwargs,
                                                                {"emodel",
                                                                "imodel"});
           }))
      ;

  py::class_<DamagedStandardKinematicModel, StandardKinematicModel,
      std::shared_ptr<DamagedStandardKinematicModel>>(m, "DamagedStandardKinematicModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
           {
            return create_object_python<DamagedStandardKinematicModel>(args, kwargs,
                                                                {"emodel",
                                                                "imodel",
                                                                "dmodel"});
           }))
      ;

} // PYBIND11_MODULE

} // namespace neml
