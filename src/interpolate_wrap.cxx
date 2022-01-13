#include "pyhelp.h" // include first to avoid annoying redef warning

#include "interpolate.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(interpolate, m) {
  py::module::import("neml.objects");

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

  py::class_<GenericPiecewiseInterpolate, Interpolate, std::shared_ptr<GenericPiecewiseInterpolate>>(m, "GenericPiecewiseInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<GenericPiecewiseInterpolate>(args, kwargs,
                                                                  {"points", "functions"});
        }))
      ;

  py::class_<PiecewiseLogLinearInterpolate, Interpolate, std::shared_ptr<PiecewiseLogLinearInterpolate>>(m, "PiecewiseLogLinearInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<PiecewiseLogLinearInterpolate>(args, kwargs,
                                                                  {"points", "values"});
        }))
      ;

  py::class_<PiecewiseSemiLogXLinearInterpolate, Interpolate,
      std::shared_ptr<PiecewiseSemiLogXLinearInterpolate>>(m, "PiecewiseSemiLogXLinearInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<PiecewiseSemiLogXLinearInterpolate>(args, kwargs,
                                                                  {"points", "values"});
        }))
      ;

  py::class_<ConstantInterpolate, Interpolate, std::shared_ptr<ConstantInterpolate>>(m, "ConstantInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ConstantInterpolate>(args, kwargs, {"v"});
        }))
      ;

  py::class_<ExpInterpolate, Interpolate, std::shared_ptr<ExpInterpolate>>(m, "ExpInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ExpInterpolate>(args, kwargs, {"A","B"});
        }))
      ;

  py::class_<PowerLawInterpolate, Interpolate,
      std::shared_ptr<PowerLawInterpolate>>(m, "PowerLawInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<PowerLawInterpolate>(args, kwargs, {"A","B"});
        }))
      ;

  py::class_<MTSShearInterpolate, Interpolate, std::shared_ptr<MTSShearInterpolate>>(m, "MTSShearInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<MTSShearInterpolate>(args, kwargs, 
                                                           {"V0", "D", "T0"});
        }))
      ;

  py::class_<MTSInterpolate, Interpolate, std::shared_ptr<MTSInterpolate>>(m, "MTSInterpolate")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<MTSInterpolate>(args, kwargs, 
                  {"tau0", "g0", "q", "p", "k", "b", "mu"});
        }))
      ;
}

} // namespace neml
