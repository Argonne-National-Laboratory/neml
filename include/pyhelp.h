#ifndef PYHELP_H
#define PYHELP_H

// To fix redef warnings
#include "Python.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "objects.h"
#include "interpolate.h"
#include "nemlerror.h"

#define PICKLEABLE(T) \
      .def(py::pickle(\
              [](std::shared_ptr<T> p) {\
                return p->serialize("object", "");\
              },\
              [](std::string state) {\
                return std::dynamic_pointer_cast<T>(get_object_string(state));\
              }\
              ))\

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

// Allocate a new, zeroed rank 3
template<class T> py::array_t<T> alloc_3d(size_t m, size_t n, size_t o)
{
  auto arr = py::array(py::buffer_info(
          nullptr,
          sizeof(T),
          py::format_descriptor<T>::value,
          3,
          {m, n, o},
          {sizeof(T) * n * o, sizeof(T) * o, sizeof(T)}
          ));
  auto ptr = arr2ptr<T>(arr);
  std::fill(ptr, ptr + m * n * o, 0);
  return arr;
}

/// Allocate an array of arbitrary size
template<class T> py::array_t<T> alloc_tensor(std::vector<size_t> shape)
{
  std::vector<size_t> stride(shape.size());
  stride[shape.size()-1] = sizeof(T);
  size_t ntotal = shape[shape.size() - 1];
  for (int i = shape.size() - 2; i >= 0; i--) {
    stride[i] = stride[i+1] * shape[i+1];
    ntotal *= shape[i];
  }

  auto arr = py::array(py::buffer_info(
          nullptr,
          sizeof(T),
          py::format_descriptor<T>::value,
          shape.size(),
          shape,
          stride));
  auto ptr = arr2ptr<T>(arr);
  std::fill(ptr, ptr + ntotal, 0);
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
      // If it's a double then we actually want a ConstantInterpolate
      // I hate using exceptions here, but what I actually want is:
      //  "is this something that can be cast to a double"
      // and not
      //  "is this actually a c++ double"
      try {
        double v = py::cast<double>(value);
        pset.assign_parameter(name, make_constant(v));
        break;
      }
      catch (py::cast_error & e) {

      }
      pset.assign_parameter(name, py::cast<std::shared_ptr<NEMLObject>>(value));
      break;
    case TYPE_VEC_NEML_OBJECT:
      // If it's a vector<double> then we actually want a vector of
      // ConstantInterpolates
      try {
        std::vector<double> v = py::cast<std::vector<double>>(value);
        std::vector<std::shared_ptr<NEMLObject>> vect;
        for (auto it = v.begin(); it != v.end(); ++it) {
          vect.push_back(make_constant(*it));
        }
        pset.assign_parameter(name, vect);
        break;
      }
      catch (py::cast_error & e) {

      }
      pset.assign_parameter(name, py::cast<std::vector<std::shared_ptr<NEMLObject>>>(value));
      break;
    case TYPE_STRING:
      pset.assign_parameter(name, py::cast<std::string>(value));
      break;
    case TYPE_SLIP:
      pset.assign_parameter(name, py::cast<list_systems>(value));
      break;
    case TYPE_TWIN:
      pset.assign_parameter(name, py::cast<twin_systems>(value));
      break;
    case TYPE_SIZE_TYPE:
      pset.assign_parameter(name, py::cast<size_t>(value));
      break;
    case TYPE_VEC_SIZE_TYPE:
      pset.assign_parameter(name, py::cast<std::vector<size_t>>(value));
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

  for (size_t i=0; i<args.size(); i++) {
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
