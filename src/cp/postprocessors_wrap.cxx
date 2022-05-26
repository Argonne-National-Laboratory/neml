#include "pyhelp.h"

#include "cp/postprocessors.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(postprocessors, m) {
  m.doc() = "Crystal plasticity postprocessors";

  py::class_<CrystalPostprocessor, NEMLObject, std::shared_ptr<CrystalPostprocessor>>(m,"CrystalPostprocessor")
      .def("populate_history", &CrystalPostprocessor::populate_hist)
      .def("init_history", &CrystalPostprocessor::init_hist)
      .def("act", &CrystalPostprocessor::act)
      ;

  py::class_<PTRTwinReorientation, CrystalPostprocessor, std::shared_ptr<PTRTwinReorientation>>(m, "PTRTwinReorientation")
    .def(py::init([](py::args args, py::kwargs kwargs)
         {
            return create_object_python<PTRTwinReorientation>(
                args, kwargs, {"threshold"});
         }))
  ;
} // PYBIND11_MODULE

} // namespace neml
