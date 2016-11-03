#include "neml.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(neml) {
  py::module m("neml", "Base class material models.");
  
  py::class_<NEMLModel, std::shared_ptr<NEMLModel>>(m, "NEMLModel")
      .def_property_readonly("nstore", &NEMLModel::nstore, "Number of variables the program needs to store.")
      .def("init_store",
           [](const NEMLModel & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nstore());
            m.init_store(arr2ptr<double>(h));
            return h;
           }, "Initialize stored variables.")
      ;

  return m.ptr();
}


} // namespace neml
