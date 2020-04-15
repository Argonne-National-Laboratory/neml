#include "pyhelp.h" // include first to avoid annoying redef warning

#include "parse.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(parse, m) {
  m.doc() = "Python wrapper to read XML input files.";
  
  m.def("parse_xml", &parse_xml);
  m.def("parse_string", &parse_string);

  py::register_exception<NodeNotFound>(m, "NodeNotFound");
  py::register_exception<DuplicateNode>(m, "DuplicateNode");
  py::register_exception<InvalidType>(m, "InvalidType");
  py::register_exception<UnknownParameterXML>(m, "UnknownParameterXML");
  py::register_exception<UnregisteredXML>(m, "UnregisteredXML");
}

} // namespace neml
