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
      .def("set_data", 
           [](History & m, py::array_t<double> arr)
           {
            m.set_data(arr2ptr<double>(arr));
           }, "Set data as a numpy array")
      .def_property_readonly("size", &History::size)
      .def_property_readonly("store", &History::store)
      .def("add_scalar", &History::add_scalar)
      .def("get_scalar", 
           [](History & m, std::string name) -> double
           {
            return m.get_scalar(name);
           }, "Get a scalar")
      .def("set_scalar",
           [](History & m, std::string name, double value)
           {
            m.get_scalar(name) = value;
           }, "Set a scalar")
      .def("add_array", &History::add_array)
      .def("array_size", &History::array_size)
      .def("get_array", 
           [](History & m, std::string name) -> py::array_t<double>
           {
            auto capsule = py::capsule(m.get_array(name), [](void *v) {;});
            return py::array_t<double>(m.array_size(name),
                                       m.get_array(name),
                                       capsule);
           }, "Get a numpy array")
      .def("add_vector", 
           [](History & m, std::string name)
           {
            m.add_object<Vector>(name);
           }, "Add a vector")
      .def("get_vector",
           [](History & m, std::string name) -> Vector
           {
            return m.get_object<Vector>(name);
           }, "Get a vector")
      .def("get_object", &History::get_object<Vector>)
      ;
}

} // namespace neml
