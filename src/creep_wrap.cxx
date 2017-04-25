#include "creep.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(creep) {
  py::module m("creep", "Creep models.");
  
  py::class_<ScalarCreepRule, std::shared_ptr<ScalarCreepRule>>(m, "ScalarCreepRule")
      .def("g",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            int ier = m.g(seq, eeq, t, T, gv);
            py_error(ier);
            return gv;
           }, "Evaluate creep rate.")

      .def("dg_ds",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            int ier = m.dg_ds(seq, eeq, t, T, gv);
            py_error(ier);
            return gv;
           }, "Evaluate creep rate derivative wrt stress.")

      .def("dg_de",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            int ier = m.dg_de(seq, eeq, t, T, gv);
            py_error(ier);
            return gv;
           }, "Evaluate creep rate derivative wrt strain.")

      .def("dg_dt",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            int ier = m.dg_dt(seq, eeq, t, T, gv);
            py_error(ier);
            return gv;
           }, "Evaluate creep rate wrt time.")

      .def("dg_dT",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            int ier = m.dg_dT(seq, eeq, t, T, gv);
            py_error(ier);
            return gv;
           }, "Evaluate creep rate wrt temperature.")
      ;
  
  py::class_<PowerLawCreep, std::shared_ptr<PowerLawCreep>>(m, "PowerLawCreep", py::base<ScalarCreepRule>())
      .def(py::init<double, double>(), py::arg("A"), py::arg("n"))
      .def(py::init<std::shared_ptr<Interpolate>, std::shared_ptr<Interpolate>>(), py::arg("A"), py::arg("n"))

      .def("A", &PowerLawCreep::A)
      .def("n", &PowerLawCreep::n)
    ;

  return m.ptr();
} // PLUGIN(creep)

} // namespace neml
