#include "../pyhelp.h" // include first to avoid redef warning

#include "singlecrystal.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(singlecrystal, m) {
  py::module::import("neml.models");

  m.doc() = "Single crystal constitutive models";

  py::class_<SingleCrystalModel, NEMLModel_ldi, std::shared_ptr<SingleCrystalModel>>(m, "SingleCrystalModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<SingleCrystalModel>(args, 
                                                                      kwargs,
                                                                      {"kinematics", "initial_rotation", "alpha"});
                    }))
    ;
}

} // namespace neml
