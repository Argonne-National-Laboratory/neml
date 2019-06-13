#include "../pyhelp.h"

#include "kinematics.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(kinematics, m) {
  m.doc() = "Kinematic models for crystal plasticity";

  py::class_<KinematicModel, std::shared_ptr<KinematicModel>>(m,
                                                              "KinematicModel")
      .def("d_p", &KinematicModel::d_p)
      .def("w_p", &KinematicModel::w_p)
      ;

  py::class_<NoInelasticity, KinematicModel, std::shared_ptr<NoInelasticity>>(m,
                                                                              "NoInelasticity")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<NoInelasticity>(args, kwargs, {});
                    }))
      ;

} // PYBIND11_MODULE

} // namespace neml
