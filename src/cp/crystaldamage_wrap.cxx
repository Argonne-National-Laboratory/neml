#include "pyhelp.h"

#include "cp/crystaldamage.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(crystaldamage, m) {
  m.doc() = "Crystal plasticity damage models";

  py::class_<CrystalDamageModel, HistoryNEMLObject,
      std::shared_ptr<CrystalDamageModel>>(m, "CrystalDamageModel")
    .def_property_readonly("nvars", &CrystalDamageModel::nvars)
    .def_property_readonly("varnames", &CrystalDamageModel::varnames)
    .def("set_varnames", &CrystalDamageModel::set_varnames)
    .def("projection", &CrystalDamageModel::projection)
    .def("d_projection_d_stress", &CrystalDamageModel::d_projection_d_stress)
    .def("d_projection_d_history", &CrystalDamageModel::d_projection_d_history)
    .def("damage_rate", &CrystalDamageModel::damage_rate)
    .def("d_damage_d_stress", &CrystalDamageModel::d_damage_d_stress)
    .def("d_damage_d_history", &CrystalDamageModel::d_damage_d_history)
  ;

  py::class_<NilDamageModel, CrystalDamageModel, 
      std::shared_ptr<NilDamageModel>>(m, "NilDamageModel")
    .def(py::init([](py::args args, py::kwargs kwargs)
                  {
                    return create_object_python<NilDamageModel>(args, 
                                                                kwargs, {});
                  }))
  ;

  py::class_<PlanarDamageModel, CrystalDamageModel, 
      std::shared_ptr<PlanarDamageModel>>(m, "PlanarDamageModel")
    .def(py::init([](py::args args, py::kwargs kwargs)
                  {
                    return create_object_python<PlanarDamageModel>(
                        args, kwargs, {"damage", "shear_transform",
                        "normal_transform", "lattice"});
                  }))
  ;

  py::class_<SlipPlaneDamage, NEMLObject,
      std::shared_ptr<SlipPlaneDamage>>(m, "SlipPlaneDamage")
    .def("setup", &SlipPlaneDamage::setup)
    .def("damage_rate", &SlipPlaneDamage::damage_rate)
    .def("d_damage_rate_d_shear", &SlipPlaneDamage::d_damage_rate_d_shear)
    .def("d_damage_rate_d_slip", &SlipPlaneDamage::d_damage_rate_d_slip)
    .def("d_damage_rate_d_normal", &SlipPlaneDamage::d_damage_rate_d_normal)
    .def("d_damage_rate_d_damage", &SlipPlaneDamage::d_damage_rate_d_damage)
  ;

  py::class_<WorkPlaneDamage, SlipPlaneDamage,
      std::shared_ptr<WorkPlaneDamage>>(m, "WorkPlaneDamage")
    .def(py::init([](py::args args, py::kwargs kwargs)
                  {
                    return create_object_python<WorkPlaneDamage>(args, kwargs,
                                                                 {});
                  }))
  ;

  py::class_<TransformationFunction, NEMLObject,
      std::shared_ptr<TransformationFunction>>(m, "TransformationFunction")
    .def("map", &TransformationFunction::map)
    .def("d_map_d_damage", &TransformationFunction::d_map_d_damage)
    .def("d_map_d_normal", &TransformationFunction::d_map_d_normal)
  ;

  py::class_<SigmoidTransformation, TransformationFunction,
      std::shared_ptr<SigmoidTransformation>>(m, "SigmoidTransformation")
    .def(py::init([](py::args args, py::kwargs kwargs)
                  {
                    return create_object_python<SigmoidTransformation>(
                        args, kwargs, {"c", "beta"});
                  }))
  ;

  py::class_<SwitchTransformation, TransformationFunction,
      std::shared_ptr<SwitchTransformation>>(m, "SwitchTransformation")
    .def(py::init([](py::args args, py::kwargs kwargs)
                  {
                    return create_object_python<SwitchTransformation>(
                        args, kwargs, {"base"});
                  }))
  ;

} // PYBIND11_MODULE

} // namespace neml
