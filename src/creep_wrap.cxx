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

  py::class_<CreepModelTrialState>(m, "CreepModelTrialState", py::base<TrialState>())
      ;

  py::class_<CreepModel, std::shared_ptr<CreepModel>>(m, "CreepModel")
      .def("update",
           [](CreepModel & m, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n) -> std::tuple<py::array_t<double>, py::array_t<double>>
           {
            auto e_np1 = alloc_vec<double>(6);
            auto A_np1 = alloc_mat<double>(6,6);

            int ier = m.update(arr2ptr<double>(s_np1), arr2ptr<double>(e_np1), arr2ptr<double>(e_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(A_np1));
            py_error(ier);

            return std::make_tuple(e_np1, A_np1);

           }, "Update to the next creep strain & tangent derivative.")

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

      
      .def("make_trial_state",
           [](CreepModel & m, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n) -> std::unique_ptr<CreepModelTrialState>
           {
            std::unique_ptr<CreepModelTrialState> ts(new CreepModelTrialState);
            int ier = m.make_trial_state(arr2ptr<double>(s_np1),
                                         arr2ptr<double>(e_n),
                                         T_np1, T_n, t_np1, t_n, *ts);
            py_error(ier);

            return ts;

           }, "Setup trial state for solve")

      // Remove if/when pybind11 supports multiple inheritance
      .def_property_readonly("nparams", &CreepModel::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](CreepModel & m, CreepModelTrialState & ts) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x), &ts);
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](CreepModel & m, py::array_t<double, py::array::c_style> x, CreepModelTrialState & ts) -> std::tuple<py::array_t<double>, py::array_t<double>>
           {
            auto R = alloc_vec<double>(m.nparams());
            auto J = alloc_mat<double>(m.nparams(), m.nparams());
            
            int ier = m.RJ(arr2ptr<double>(x), &ts, arr2ptr<double>(R), arr2ptr<double>(J));
            py_error(ier);

            return std::make_tuple(R, J);
           }, "Residual and jacobian.")
      ;
      // End remove block


    ;

  py::class_<J2CreepModel, std::shared_ptr<J2CreepModel>>(m, "J2CreepModel", py::base<CreepModel>())
      .def(py::init<std::shared_ptr<ScalarCreepRule>, double, int , bool>(),
           py::arg("rule"), py::arg("tol") = 1.0e-10, py::arg("miter") = 25,
           py::arg("verbose") = false)
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

  py::class_<RegionKMCreep, std::shared_ptr<RegionKMCreep>>(m, "RegionKMCreep", py::base<ScalarCreepRule>())
      .def(py::init<std::vector<double>, std::vector<double>, std::vector<double>, double, double, double, std::shared_ptr<LinearElasticModel>>(), py::arg("cuts"), py::arg("A"), py::arg("B"), py::arg("kboltz"), py::arg("b"), py::arg("eps0"), py::arg("emodel"))
  ;

  py::class_<NortonBaileyCreep, std::shared_ptr<NortonBaileyCreep>>(m, "NortonBaileyCreep", py::base<ScalarCreepRule>())
      .def(py::init<double, double, double>(), py::arg("A"), py::arg("m"), py::arg("n"))
      .def(py::init<std::shared_ptr<Interpolate>, std::shared_ptr<Interpolate>, std::shared_ptr<Interpolate>>(), py::arg("A"), py::arg("m"), py::arg("n"))

      .def("A", &NortonBaileyCreep::A)
      .def("m", &NortonBaileyCreep::m)
      .def("n", &NortonBaileyCreep::n)
    ;

  py::class_<MukherjeeCreep, std::shared_ptr<MukherjeeCreep>>(m, "MukherjeeCreep", py::base<ScalarCreepRule>())
      .def(py::init<std::shared_ptr<LinearElasticModel>, double, double, double, double, double, double, double>(), 
           py::arg("emodel"), py::arg("A"), py::arg("n"), py::arg("D0"), py::arg("Q"), py::arg("b"),
           py::arg("k"), py::arg("R"))

      .def_property_readonly("A", &MukherjeeCreep::A)
      .def_property_readonly("n", &MukherjeeCreep::n)
      .def_property_readonly("D0", &MukherjeeCreep::D0)
      .def_property_readonly("Q", &MukherjeeCreep::Q)
      .def_property_readonly("b", &MukherjeeCreep::b)
      .def_property_readonly("k", &MukherjeeCreep::k)
      .def_property_readonly("R", &MukherjeeCreep::R)
    ;

  return m.ptr();
} // PLUGIN(creep)

} // namespace neml
