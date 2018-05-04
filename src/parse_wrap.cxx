#include "parse.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_MODULE(parse, m) {
  m.doc() = "Python wrapper to read XML input files.";
  
  m.def("parse_xml", 
        [](std::string fname, std::string mname) -> std::shared_ptr<NEMLModel>
        {
          int ier;
          std::unique_ptr<NEMLModel> m = parse_xml(fname, mname, ier);
          py_error(ier);

          return std::move(m);
        }, "Return a model from an XML description");
}

} // namespace neml
