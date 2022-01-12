#include "pyhelp.h" // include first to avoid annoying redef warning

#include "objects.h"

// #include "parse.h"
// #include "deparse.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(objects, m) {
  m.doc() = "Factory system for creating objects.";

  py::class_<NEMLObject, std::shared_ptr<NEMLObject>>(m, "NEMLObject")
      .def("serialize", &NEMLObject::serialize, py::arg("top_name") = "object",
           py::arg("top_node") = "")
      // This should work but doesn't.  Instead we have the stupid
      // PICKLEABLE(T) macro in pyhelp.h at the fully-defined class level
      /*
      .def(py::pickle(
              [](std::shared_ptr<NEMLObject> p) {
                return p->serialize("object", "");
              },
              [](std::string state) {
                return get_object_string(state);
              }
              ))
      */
      ;
}

} // namespace neml
