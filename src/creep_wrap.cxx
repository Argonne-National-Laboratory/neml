#include "pyhelp.h" // include first to avoid annoying redef warning

#include "creep.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(creep, m) {
  py::module::import("neml.objects");
  py::module::import("neml.solvers");

  m.doc() = "Separate creep models to combine with base NEML models";

  py::class_<CreepModelTrialState, TrialState>(m, "CreepModelTrialState")
      ;

  py::class_<CreepModel, NEMLObject, Solvable, std::shared_ptr<CreepModel>>(m, "CreepModel")
      .def("update",
           [](CreepModel & m, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n) -> std::tuple<py::array_t<double>, py::array_t<double>>
           {
            auto e_np1 = alloc_vec<double>(6);
            auto A_np1 = alloc_mat<double>(6,6);

            m.update(arr2ptr<double>(s_np1), arr2ptr<double>(e_np1), arr2ptr<double>(e_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(A_np1));

            return std::make_tuple(e_np1, A_np1);

           }, "Update to the next creep strain & tangent derivative.")

      .def("f",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto fv = alloc_vec<double>(6);
            m.f(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(fv));
            return fv;
           }, "Evaluate creep rate.")

      .def("df_ds",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto dfv = alloc_mat<double>(6,6);
            m.df_ds(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(dfv));
            return dfv;
           }, "Evaluate creep rate derivative wrt stress.")

      .def("df_de",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto dfv = alloc_mat<double>(6,6);
            m.df_de(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(dfv));
            return dfv;
           }, "Evaluate creep rate derivative wrt strain.")

      .def("df_dt",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto dfv = alloc_vec<double>(6);
            m.df_dt(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(dfv));
            return dfv;
           }, "Evaluate creep rate derivative wrt time.")

      .def("df_dT",
           [](const CreepModel & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> e, double t, double T) -> py::array_t<double>
           {
            auto dfv = alloc_vec<double>(6);
            m.df_dT(arr2ptr<double>(s), arr2ptr<double>(e), t, T, arr2ptr<double>(dfv));
            return dfv;
           }, "Evaluate creep rate derivative wrt temperature.")

      
      .def("make_trial_state",
           [](CreepModel & m, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n) -> std::unique_ptr<CreepModelTrialState>
           {
            std::unique_ptr<CreepModelTrialState> ts(new CreepModelTrialState);
            m.make_trial_state(arr2ptr<double>(s_np1),
                                         arr2ptr<double>(e_n),
                                         T_np1, T_n, t_np1, t_n, *ts);

            return ts;

           }, "Setup trial state for solve")
    ;

  py::class_<J2CreepModel, CreepModel, std::shared_ptr<J2CreepModel>>(m, "J2CreepModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<J2CreepModel>(args, kwargs, {"rule"});
        }))
    ;

  py::class_<ScalarCreepRule, NEMLObject, std::shared_ptr<ScalarCreepRule>>(m, "ScalarCreepRule")
      .def("g",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            m.g(seq, eeq, t, T, gv);
            return gv;
           }, "Evaluate creep rate.")

      .def("dg_ds",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            m.dg_ds(seq, eeq, t, T, gv);
            return gv;
           }, "Evaluate creep rate derivative wrt stress.")

      .def("dg_de",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            m.dg_de(seq, eeq, t, T, gv);
            return gv;
           }, "Evaluate creep rate derivative wrt strain.")

      .def("dg_dt",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            m.dg_dt(seq, eeq, t, T, gv);
            return gv;
           }, "Evaluate creep rate wrt time.")

      .def("dg_dT",
           [](const ScalarCreepRule & m, double seq, double eeq, double t, double T) -> double
           {
            double gv;
            m.dg_dT(seq, eeq, t, T, gv);
            return gv;
           }, "Evaluate creep rate wrt temperature.")
      ;
  
  py::class_<PowerLawCreep, ScalarCreepRule, std::shared_ptr<PowerLawCreep>>(m, "PowerLawCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<PowerLawCreep>(args, kwargs, {"A", "n"});
        }))

      .def("A", &PowerLawCreep::A)
      .def("n", &PowerLawCreep::n)
    ;

  py::class_<NormalizedPowerLawCreep, ScalarCreepRule, std::shared_ptr<NormalizedPowerLawCreep>>(m, "NormalizedPowerLawCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<NormalizedPowerLawCreep>(args, kwargs, {"s0", "n"});
        }))
    ;

  py::class_<BlackburnMinimumCreep, ScalarCreepRule, std::shared_ptr<BlackburnMinimumCreep>>(m, "BlackburnMinimumCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<BlackburnMinimumCreep>(args, kwargs, {"A","n","beta","R","Q"});
        }))
    ;

  py::class_<SwindemanMinimumCreep, ScalarCreepRule, std::shared_ptr<SwindemanMinimumCreep>>(m, "SwindemanMinimumCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<SwindemanMinimumCreep>(args, kwargs, {"C", "n", "V", "Q"});
        }))
    ;

  py::class_<RegionKMCreep, ScalarCreepRule, std::shared_ptr<RegionKMCreep>>(m, "RegionKMCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<RegionKMCreep>(args, kwargs, {"cuts", "A", "B",
                                                     "kboltz", "b", "eps0", "emodel"});
        }))
  ;

  py::class_<NortonBaileyCreep, ScalarCreepRule, std::shared_ptr<NortonBaileyCreep>>(m, "NortonBaileyCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<NortonBaileyCreep>(args, kwargs, {"A", "m",
                                                         "n"});
        }))
      .def("A", &NortonBaileyCreep::A)
      .def("m", &NortonBaileyCreep::m)
      .def("n", &NortonBaileyCreep::n)
    ;

  py::class_<MukherjeeCreep, ScalarCreepRule, std::shared_ptr<MukherjeeCreep>>(m, "MukherjeeCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<MukherjeeCreep>(args, kwargs, {"emodel", "A",
                                                      "n", "D0", "Q", "b", "k", "R"});
        }))
      .def_property_readonly("A", &MukherjeeCreep::A)
      .def_property_readonly("n", &MukherjeeCreep::n)
      .def_property_readonly("D0", &MukherjeeCreep::D0)
      .def_property_readonly("Q", &MukherjeeCreep::Q)
      .def_property_readonly("b", &MukherjeeCreep::b)
      .def_property_readonly("k", &MukherjeeCreep::k)
      .def_property_readonly("R", &MukherjeeCreep::R)
    ;

  py::class_<GenericCreep, ScalarCreepRule, std::shared_ptr<GenericCreep>>(m, "GenericCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<GenericCreep>(args, kwargs, {"cfn"});
        }))
      ;

  py::class_<MinCreep225Cr1MoCreep, ScalarCreepRule, std::shared_ptr<MinCreep225Cr1MoCreep>>(m, "MinCreep225Cr1MoCreep")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<MinCreep225Cr1MoCreep>(args, kwargs, {});
                    }))
  ;
} // MODULE(creep, m)

} // namespace neml
