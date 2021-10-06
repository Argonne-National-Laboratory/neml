#include "pyhelp.h"

#include "math/projections.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(projections, m) {
  m.doc() = "Projection operators onto planes";

  m.def("normal_projection", &normal_projection);
  m.def("normal_projection_ss", &normal_projection_ss);

  m.def("shear_projection", &shear_projection);
  m.def("shear_projection_ss", &shear_projection_ss);

} // MODULE end

} // namespace neml
