#include "../pyhelp.h"

#include "crystaldamage.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(crystaldamage, m) {
  m.doc() = "Crystal plasticity damage models";

  py::class_<CrystalDamageModel, NEMLObject,
      std::shared_ptr<CrystalDamageModel>>(m, "CrystalDamageModel")
    .def_property_readonly("nvars", &CrystalDamageModel::nvars)
    .def_property_readonly("varnames", &CrystalDamageModel::varnames)
    .def("set_varnames", &CrystalDamageModel::set_varnames)
    .def("setup", &CrystalDamageModel::setup)
    .def("populate_history", &CrystalDamageModel::populate_history)
    .def("init_history", &CrystalDamageModel::init_history)
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

} // PYBIND11_MODULE

} // namespace neml
