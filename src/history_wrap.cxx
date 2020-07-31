#include "pyhelp.h" // include first to avoid annoying redef warning

#include "history.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(history, m) {
  m.doc() = "Internal variable tracking system.";

  py::class_<History, std::shared_ptr<History>>(m, "History",
                                                py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<bool>(), py::arg("store"))
      .def_buffer(
          [](History & m) -> py::buffer_info 
          {
            return py::buffer_info(
                m.rawptr(),
                sizeof(double),
                py::format_descriptor<double>::format(),
                1,
                {m.size()},
                {sizeof(double)});
          })
      .def("deepcopy", &History::deepcopy)
      .def("set_data", 
           [](History & m, py::array_t<double> arr)
           {
            m.set_data(arr2ptr<double>(arr));
           }, "Set data as a numpy array")
      .def("copy_data",
           [](History & m, py::array_t<double> arr)
           {
            m.copy_data(arr2ptr<double>(arr));
           }, "Copy data in from an array")
      .def_property_readonly("size", &History::size)
      .def_property_readonly("store", &History::store)
      .def_property_readonly("items", &History::items)
      .def("add_scalar", 
           [](History & m, std::string name)
              {
                m.add<double>(name);
              }, "Add a scalar")
      .def("get_scalar", 
           [](History & m, std::string name) -> double
           {
            return m.get<double>(name);
           }, "Get a scalar")
      .def("set_scalar",
           [](History & m, std::string name, double value)
           {
            m.get<double>(name) = value;
           }, "Set a scalar")
      .def("add_vector", 
           [](History & m, std::string name)
           {
            m.add<Vector>(name);
           }, "Add a vector")
      .def("get_vector",
           [](History & m, std::string name) -> Vector
           {
            return m.get<Vector>(name);
           }, "Get a vector")
      .def("set_vector",
           [](History & m, std::string name, Vector v)
           {
            m.get<Vector>(name) = v;
           }, "Set a vector")
      .def("add_ranktwo", 
           [](History & m, std::string name)
           {
            m.add<RankTwo>(name);
           }, "Add a vector")
      .def("get_ranktwo",
           [](History & m, std::string name) -> RankTwo
           {
            return m.get<RankTwo>(name);
           }, "Get a general tensor")
      .def("set_ranktwo",
           [](History & m, std::string name, RankTwo v)
           {
            m.get<RankTwo>(name) = v;
           }, "Set a general tensor")

      .def("add_symmetric", 
           [](History & m, std::string name)
           {
            m.add<Symmetric>(name);
           }, "Add a vector")
      .def("get_symmetric",
           [](History & m, std::string name) -> Symmetric
           {
            return m.get<Symmetric>(name);
           }, "Get a symmetric tensor")
      .def("set_symmetric",
           [](History & m, std::string name, Symmetric v)
           {
            m.get<Symmetric>(name) = v;
           }, "Set a symmetric tensor")
      .def("add_skew", 
           [](History & m, std::string name)
           {
            m.add<Skew>(name);
           }, "Add a skew")
      .def("get_skew",
           [](History & m, std::string name) -> Skew
           {
            return m.get<Skew>(name);
           }, "Get a skew tensor")
      .def("set_skew",
           [](History & m, std::string name, Skew v)
           {
            m.get<Skew>(name) = v;
           }, "Set a skew")
      .def("add_orientation", 
           [](History & m, std::string name)
           {
            m.add<Orientation>(name);
           }, "Add an orientation")
      .def("get_orientation",
           [](History & m, std::string name) -> Orientation
           {
            return m.get<Orientation>(name);
           }, "Get an orientation")
      .def("set_orientation",
           [](History & m, std::string name, Orientation v)
           {
            m.get<Orientation>(name) = v;
           }, "Set an orientation")
        .def("scalar_multiply", &History::scalar_multiply)
        .def(py::self += py::self)
        .def("zero", &History::zero)
        .def("split", &History::split, py::arg("group"), py::arg("after") = true)
        .def("add_union", &History::add_union)
        .def("contains", &History::contains)
        .def("subset", &History::subset)
        .def("reorder", &History::reorder)
      ;
}

} // namespace neml
