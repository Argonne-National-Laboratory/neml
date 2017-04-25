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
 
  py::class_<CreepModel, std::shared_ptr<CreepModel>>(m, "CreepModel")
      .def("f",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto fv = alloc_vec<double>(6);
            int ier = m.f(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(fv));
            py_error(ier);
            return fv;
           }, "Evaluate creep rate.")

      .def("df_ds",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto dfv = alloc_mat<double>(6,6);
            int ier = m.df_ds(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(dfv));
            py_error(ier);
            return dfv;
           }, "Evaluate creep rate derivative wrt stress.")

      .def("df_de",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto dfv = alloc_mat<double>(6,6);
            int ier = m.df_de(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(dfv));
            py_error(ier);
            return dfv;
           }, "Evaluate creep rate derivative wrt strain.")

      .def("df_dt",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto dfv = alloc_vec<double>(6);
            int ier = m.df_dt(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(dfv));
            py_error(ier);
            return dfv;
           }, "Evaluate creep rate derivative wrt time.")

      .def("df_dT",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto dfv = alloc_vec<double>(6);
            int ier = m.df_dT(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(dfv));
            py_error(ier);
            return dfv;
           }, "Evaluate creep rate derivative wrt temperature.")
    ;

  py::class_<J2CreepModel, std::shared_ptr<J2CreepModel>>(m, "J2CreepModel", py::base<CreepModel>())
      .def(py::init<std::shared_ptr<ScalarCreepRule>>(), py::arg("rule"))
    ;

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

  py::class_<NortonBaileyCreep, std::shared_ptr<NortonBaileyCreep>>(m, "NortonBaileyCreep", py::base<ScalarCreepRule>())
      .def(py::init<double, double, double>(), py::arg("A"), py::arg("m"), py::arg("n"))
      .def(py::init<std::shared_ptr<Interpolate>, std::shared_ptr<Interpolate>, std::shared_ptr<Interpolate>>(), py::arg("A"), py::arg("m"), py::arg("n"))

      .def("A", &NortonBaileyCreep::A)
      .def("m", &NortonBaileyCreep::m)
      .def("n", &NortonBaileyCreep::n)
    ;

  return m.ptr();
} // PLUGIN(creep)

} // namespace neml
