#include "pyhelp.h"

#include "interpolate.h"

namespace neml {

void assign_python_parameter(ParameterSet & pset, std::string name,
                             py::object value)
{
  switch (pset.get_object_type(name)) {
    case TYPE_DOUBLE:
      pset.assign_parameter(name, py::cast<double>(value));
      break;
    case TYPE_INT:
      pset.assign_parameter(name, py::cast<int>(value));
      break;
    case TYPE_BOOL:
      pset.assign_parameter(name, py::cast<bool>(value));
      break;
    case TYPE_VEC_DOUBLE:
      pset.assign_parameter(name, py::cast<std::vector<double>>(value));
      break;
    case TYPE_NEML_OBJECT:
      pset.assign_parameter(name, py::cast<std::shared_ptr<NEMLObject>>(value));
      break;
    case TYPE_VEC_NEML_OBJECT:
      pset.assign_parameter(name, py::cast<std::vector<std::shared_ptr<NEMLObject>>>(value));
      break;
    case TYPE_INTERPOLATE:
      try {
        pset.assign_parameter(name, py::cast<std::shared_ptr<Interpolate>>(value));
      }
      catch (py::cast_error e) {
        try {
          double v = py::cast<double>(value);
          pset.assign_parameter(name, std::make_shared<ConstantInterpolate>(v));
        }
        catch (py::cast_error e) {
          throw std::runtime_error("Cannot convert object to an Interpolate object!");
        }
      }
      break;
    case TYPE_VEC_INTERPOLATE:
      try {
        pset.assign_parameter(name, py::cast<std::vector<std::shared_ptr<Interpolate>>>(value));
      }
      catch (py::cast_error e) {
        try {
          std::vector<double> v = py::cast<std::vector<double>>(value);
          std::vector<std::shared_ptr<Interpolate>> vect;

          pset.assign_parameter(name, vect);
        }
        catch (py::cast_error e) {
          throw std::runtime_error("Cannot convert object to a vector of Interpolate objects!");
        }
      }
      break;
    default:
      throw std::runtime_error("Unrecognized object type!");
      break;
  }
}


} // namespace neml
