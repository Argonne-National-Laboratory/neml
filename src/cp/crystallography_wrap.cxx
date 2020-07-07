#include "../pyhelp.h"

#include "crystallography.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(crystallography, m) {
  py::module::import("neml.objects");

  m.doc() = "Various routines for crystallography";
  
  m.def("symmetry_rotations", &symmetry_rotations);

  py::class_<SymmetryGroup, NEMLObject, std::shared_ptr<SymmetryGroup>>(m, "SymmetryGroup")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<SymmetryGroup>(
                          args, kwargs, {"sclass"});
                    }))
      .def_property_readonly("ops", &SymmetryGroup::ops)
      .def("misorientation", &SymmetryGroup::misorientation)
      .def("misorientation_block", &SymmetryGroup::misorientation_block)
      ;

  py::class_<Lattice, NEMLObject, std::shared_ptr<Lattice>>(m, "Lattice")
      .def(py::init<Vector, Vector, Vector, std::shared_ptr<SymmetryGroup>>(),
           py::arg("a1"), py::arg("a2"), py::arg("a3"), py::arg("symmetry"))
      .def_property_readonly("a1", &Lattice::a1)
      .def_property_readonly("a2", &Lattice::a2)
      .def_property_readonly("a3", &Lattice::a3)
      .def_property_readonly("b1", &Lattice::b1)
      .def_property_readonly("b2", &Lattice::b2)
      .def_property_readonly("b3", &Lattice::b3)

      .def_property_readonly("burgers_vectors", &Lattice::burgers_vectors)
      .def_property_readonly("slip_directions", &Lattice::slip_directions)
      .def_property_readonly("slip_planes", &Lattice::slip_planes)

      .def_property_readonly("symmetry", &Lattice::symmetry)

      .def("miller2cart_direction", &Lattice::miller2cart_direction)
      .def("miller2cart_plane", &Lattice::miller2cart_plane)
      .def("equivalent_vectors", &Lattice::equivalent_vectors)
      .def("equivalent_vectors_bidirectional", &Lattice::equivalent_vectors_bidirectional)

      .def("add_slip_system", &Lattice::add_slip_system)
      .def_property_readonly("ntotal", &Lattice::ntotal)
      .def_property_readonly("ngroup", &Lattice::ngroup)
      .def("nslip", &Lattice::nslip)
      .def("M", &Lattice::M)
      .def("N", &Lattice::N)
      .def("shear", &Lattice::shear)
      .def("d_shear", &Lattice::d_shear)
      ;

  py::class_<CubicLattice, Lattice, std::shared_ptr<CubicLattice>>(m, "CubicLattice")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<CubicLattice>(
                          args, kwargs, {"a"});
                    }))
      ;

} // PYBIND11_MODULE(crystallography, m)

} // namespace neml
