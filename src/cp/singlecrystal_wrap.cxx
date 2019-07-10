#include "../pyhelp.h" // include first to avoid redef warning

#include "singlecrystal.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(singlecrystal, m) {
  py::module::import("neml.models");

  m.doc() = "Single crystal constitutive models";

  py::class_<SCTrialState, TrialState>(m, "SCTrialState")
      .def(py::init<Symmetric&,Skew&,Symmetric&,History&,Orientation&,Lattice&,double,double>())
      ;

  py::class_<SingleCrystalModel, NEMLModel_ldi, Solvable, std::shared_ptr<SingleCrystalModel>>(m, "SingleCrystalModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<SingleCrystalModel>(args, 
                                                                      kwargs,
                                                                      {"kinematics", "lattice"});
                    }))
      .def("populate_history", &SingleCrystalModel::populate_history)
      .def("init_history", &SingleCrystalModel::init_history)
      .def("get_passive_orientation", 
           [](SingleCrystalModel & m, const History & hist) -> Orientation
           {
            return m.get_passive_orientation(hist);
           }, "Get the orientation as an passive rotation (sample -> crystal)")
      .def("get_passive_orientation", 
           [](SingleCrystalModel & m, py::array_t<double, py::array::c_style>
              hist) -> Orientation
           {
            return m.get_passive_orientation(arr2ptr<double>(hist));
           }, "Get the orientation as an passive rotation (sample -> crystal)")

      .def("get_active_orientation", 
           [](SingleCrystalModel & m, const History & hist) -> Orientation
           {
            return m.get_active_orientation(hist);
           }, "Get the orientation as an active rotation (crystal -> sample)")
      .def("get_active_orientation", 
           [](SingleCrystalModel & m, py::array_t<double, py::array::c_style>
              hist) -> Orientation
           {
            return m.get_active_orientation(arr2ptr<double>(hist));
           }, "Get the orientation as an active rotation (crystal -> sample)")
      .def("set_passive_orientation",
           [](SingleCrystalModel & m, History & hist, const Orientation & q)
           {
            m.set_passive_orientation(hist, q);
           }, "Set the orientation using a passive rotation (sample -> crystal)")
      .def("set_passive_orientation",
           [](SingleCrystalModel & m, py::array_t<double, py::array::c_style>
              hist, const Orientation & q)
           {
            m.set_passive_orientation(arr2ptr<double>(hist), q);
           }, "Set the orientation using a passive rotation (sample -> crystal)")
      .def("set_active_orientation",
           [](SingleCrystalModel & m, History & hist, const Orientation & q)
           {
            m.set_active_orientation(hist, q);
           }, "Set the orientation using a active rotation (crystal -> sample)")
      .def("set_active_orientation",
           [](SingleCrystalModel & m, py::array_t<double, py::array::c_style>
              hist, const Orientation & q)
           {
            m.set_active_orientation(arr2ptr<double>(hist), q);
           }, "Set the orientation using a active rotation (crystal -> sample)")
    ;
}

} // namespace neml
