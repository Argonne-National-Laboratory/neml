#include "pyhelp.h" // include first to avoid redef warning

#include "cp/polycrystal.h"

#include "parse.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(polycrystal, m) {
  py::module::import("neml.cp.singlecrystal");

  m.doc() = "Polycrystal constitutive models";
  
  py::class_<PolycrystalModel, NEMLModel_ldi, std::shared_ptr<PolycrystalModel>>(m, "PolycrystalModel")
      .def_property_readonly("n", &PolycrystalModel::n)
      .def("orientations", 
           [](PolycrystalModel & m, py::array_t<double, py::array::c_style> h) -> std::vector<Orientation>
           {
            return m.orientations(arr2ptr<double>(h));
           }, "Return the current vector of orientations in the passive convention")
      .def("orientations_active", 
           [](PolycrystalModel & m, py::array_t<double, py::array::c_style> h) -> std::vector<Orientation>
           {
            return m.orientations_active(arr2ptr<double>(h));
           }, "Return the current vector of orientations in the active convention")
    ;

  py::class_<TaylorModel, PolycrystalModel, std::shared_ptr<TaylorModel>>(m, "TaylorModel")
      PICKLEABLE(TaylorModel)
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<TaylorModel>(args,
                                                               kwargs,
                                                               {"model", "qs"});
                    }))
    ;
}

}
