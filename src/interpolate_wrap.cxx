#include "pyhelp.h" // include first to avoid annoying redef warning

#include "interpolate.h"

#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_MODULE(interpolate, m) {
  m.doc() = "Interpolation schemes used to define model parameters.";
  
  py::class_<Interpolate, NEMLObject, std::shared_ptr<Interpolate>>(m, "Interpolate")
      .def("value", &Interpolate::value, "Interpolate to x")
      .def("derivative", &Interpolate::derivative, "Derivative at x")
      .def("__call__", 
           [](Interpolate & m, double x) -> double
           {
            return m(x);
           }, "Operator overload ()")
      .def_property_readonly("valid", &Interpolate::valid)
      ;

  py::class_<InvalidInterpolate, Interpolate, std::shared_ptr<InvalidInterpolate>>(m, "InvalidInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<InvalidInterpolate>(args, kwargs, {});
        }))
      ;

  py::class_<PolynomialInterpolate, Interpolate, std::shared_ptr<PolynomialInterpolate>>(m, "PolynomialInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<PolynomialInterpolate>(args, kwargs, 
                                                             {"coefs"});
        }))
      ;

  py::class_<PiecewiseLinearInterpolate, Interpolate, std::shared_ptr<PiecewiseLinearInterpolate>>(m, "PiecewiseLinearInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<PiecewiseLinearInterpolate>(args, kwargs,
                                                                  {"points", "values"});
        }))
      ;

  py::class_<PiecewiseLogLinearInterpolate, Interpolate, std::shared_ptr<PiecewiseLogLinearInterpolate>>(m, "PiecewiseLogLinearInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<PiecewiseLogLinearInterpolate>(args, kwargs,
                                                                  {"points", "values"});
        }))
      ;

  py::class_<ConstantInterpolate, Interpolate, std::shared_ptr<ConstantInterpolate>>(m, "ConstantInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ConstantInterpolate>(args, kwargs, {"v"});
        }))
      ;

  py::class_<MTSShearInterpolate, Interpolate, std::shared_ptr<MTSShearInterpolate>>(m, "MTSShearInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<MTSShearInterpolate>(args, kwargs, 
                                                           {"V0", "D", "T0"});
        }))
      ;
}

} // namespace neml
