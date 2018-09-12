#ifndef PYHELP_H
#define PYHELP_H

// To fix redef warnings
#include "Python.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

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

/// Map a python object into a parameter set
void assign_python_parameter(ParameterSet & pset, std::string name, 
                             py::object value);

} // namespace neml

#endif // namespace neml
