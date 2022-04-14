#include "pyhelp.h" // include first to avoid annoying redef warning

#include "parse.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(parse, m) {
  m.doc() = "Python wrapper to read XML input files.";
  
  m.def("parse_xml", &parse_xml);
  m.def("parse_string", &parse_string);

  m.def("get_object_string", &get_object_string);

  py::register_exception<InvalidType>(m, "InvalidType");
  py::register_exception<UnknownParameterXML>(m, "UnknownParameterXML");
  py::register_exception<UnregisteredXML>(m, "UnregisteredXML");
  py::register_exception<ModelNotFound>(m, "ModelNotFound");
}

} // namespace neml
