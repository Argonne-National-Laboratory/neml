#include "pyhelp.h" // include first to avoid annoying redef warning

#include "objects.h"

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(objects, m) {
  m.doc() = "Factory system for creating objects.";

  py::class_<NEMLObject, std::shared_ptr<NEMLObject>>(m, "NEMLObject")
      ;
}

} // namespace neml
