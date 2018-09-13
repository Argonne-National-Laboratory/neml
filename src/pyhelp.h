#ifndef PYHELP_H
#define PYHELP_H

// To fix redef warnings
#include "Python.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "objects.h"

namespace py = pybind11;

namespace neml {

// Convert an array to a pointer
template<class T> T* arr2ptr(py::array_t<T, py::array::c_style> arr)
{
  return static_cast<T*>(arr.request().ptr);
}

// Allocate a new, zeroed vector
template<class T> py::array_t<T> alloc_vec(size_t n)
{
  auto arr = py::array(py::buffer_info(
          nullptr,
          sizeof(T),
          py::format_descriptor<T>::value,
          1,
          {n},
          {sizeof(T)}
          ));
  auto ptr = arr2ptr<T>(arr);
  std::fill(ptr, ptr + n, 0);
  return arr;
}

// Allocate a new, zeroed matrix
template<class T> py::array_t<T> alloc_mat(size_t m, size_t n)
{
  auto arr = py::array(py::buffer_info(
          nullptr,
          sizeof(T),
          py::format_descriptor<T>::value,
          2,
          {m, n},
          {sizeof(T) * n, sizeof(T)}
          ));
  auto ptr = arr2ptr<T>(arr);
  std::fill(ptr, ptr + m * n, 0);
  return arr;
}

/// Map a python object into a parameter from a set
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

/// Create an object from args and kwargs
template<typename T>
std::shared_ptr<T> create_object_python(py::args args, py::kwargs kwargs,
                                        std::vector<std::string> names)
{
  ParameterSet pset = Factory::Creator()->provide_parameters(T::type());

  // The parameter names must map to each required arg
  if (args.size() != names.size()) {
    throw std::runtime_error("Each arg in args does not have a name in names.");
  }

  for (int i=0; i<args.size(); i++) {
    assign_python_parameter(pset, names[i], args[i]);
  }

  for (auto item : kwargs) {
    assign_python_parameter(pset, item.first.cast<std::string>(),
                            item.second.cast<py::object>());
  }

  return Factory::Creator()->create<T>(pset);
}

} // namespace neml

#endif // namespace neml
