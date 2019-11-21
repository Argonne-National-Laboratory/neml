#include "pyhelp.h" // include first to avoid annoying redef warning

#include "damage.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(damage, m) {
  py::module::import("neml.models");
  py::module::import("neml.objects");
  py::module::import("neml.solvers");
  py::module::import("neml.larsonmiller");

  m.doc() = "NEML damage models.";

  py::class_<NEMLDamagedModel_sd, NEMLModel_sd, std::shared_ptr<NEMLDamagedModel_sd>>(m, "NEMLDamagedModel_sd")
      .def_property_readonly("ndamage", &NEMLDamagedModel_sd::ndamage, "Number of damage variables.")
      .def("init_damage",
           [](NEMLDamagedModel_sd & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.ndamage());
            int ier = m.init_damage(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize damage variables.")
      ;

  py::class_<SDTrialState, TrialState>(m, "SDTrialState")
      ;

  py::class_<NEMLScalarDamagedModel_sd, NEMLDamagedModel_sd, Solvable, std::shared_ptr<NEMLScalarDamagedModel_sd>>(m, "NEMLScalarDamagedModel_sd")
      .def("damage",
           [](NEMLScalarDamagedModel_sd & m, double d_np1, double d_n, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> s_n, double T_np1, double T_n, double t_np1, double t_n) -> double
           {
            double damage;
            int ier = m.damage(d_np1, d_n, arr2ptr<double>(e_np1), arr2ptr<double>(e_n),
                     arr2ptr<double>(s_np1), arr2ptr<double>(s_n),
                     T_np1, T_n,
                     t_np1, t_n, &damage);
            py_error(ier);
            return damage;
           }, "The damage evolution equation.")
      .def("ddamage_dd",
           [](NEMLScalarDamagedModel_sd & m, double d_np1, double d_n, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> s_n, double T_np1, double T_n, double t_np1, double t_n) -> double
           {
            double ddamage;
            int ier = m.ddamage_dd(d_np1, d_n, arr2ptr<double>(e_np1), arr2ptr<double>(e_n),
                     arr2ptr<double>(s_np1), arr2ptr<double>(s_n),
                     T_np1, T_n,
                     t_np1, t_n, &ddamage);
            py_error(ier);
            return ddamage;
           }, "The derivative of the damage evolution equation wrt. damage.")
      .def("ddamage_de",
           [](NEMLScalarDamagedModel_sd & m, double d_np1, double d_n, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> s_n, double T_np1, double T_n, double t_np1, double t_n) -> py::array_t<double>
           {
            auto ddamage = alloc_vec<double>(6);
            int ier = m.ddamage_de(d_np1, d_n, arr2ptr<double>(e_np1), arr2ptr<double>(e_n),
                     arr2ptr<double>(s_np1), arr2ptr<double>(s_n),
                     T_np1, T_n,
                     t_np1, t_n, arr2ptr<double>(ddamage));
            py_error(ier);
            return ddamage;
           }, "The derivative of the damage evolution equation wrt. strain.") 
      .def("ddamage_ds",
           [](NEMLScalarDamagedModel_sd & m, double d_np1, double d_n, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> s_n, double T_np1, double T_n, double t_np1, double t_n) -> py::array_t<double>
           {
            auto ddamage = alloc_vec<double>(6);
            int ier = m.ddamage_ds(d_np1, d_n, arr2ptr<double>(e_np1), arr2ptr<double>(e_n),
                     arr2ptr<double>(s_np1), arr2ptr<double>(s_n),
                     T_np1, T_n,
                     t_np1, t_n, arr2ptr<double>(ddamage));
            py_error(ier);
            return ddamage;
           }, "The derivative of the damage evolution equation wrt. stress.") 
      .def("make_trial_state",
           [](NEMLScalarDamagedModel_sd & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n, double u_n, double p_n) -> std::unique_ptr<SDTrialState>
           {
            std::unique_ptr<SDTrialState> ts(new SDTrialState);
            int ier = m.make_trial_state(arr2ptr<double>(e_np1),
                                         arr2ptr<double>(e_n),
                                         T_np1, T_n, t_np1, t_n,
                                         arr2ptr<double>(s_n),
                                         arr2ptr<double>(h_n),
                                         u_n, p_n, *ts);
            py_error(ier);

            return ts;
           }, "Make a trial state, mostly for testing.")
      ;

  py::class_<CombinedDamageModel_sd, NEMLScalarDamagedModel_sd, std::shared_ptr<CombinedDamageModel_sd>>(m, "CombinedDamageModel_sd")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<CombinedDamageModel_sd>(args, kwargs,
                                                              {"elastic",
                                                              "models",
                                                              "base"});
        }))
      ;

  py::class_<ClassicalCreepDamageModel_sd, NEMLScalarDamagedModel_sd, std::shared_ptr<ClassicalCreepDamageModel_sd>>(m, "ClassicalCreepDamageModel_sd")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ClassicalCreepDamageModel_sd>(args, kwargs,
                                                                    {"elastic",
                                                                    "A", "xi",
                                                                    "phi", 
                                                                    "base"});
        }))
      ;

  py::class_<EffectiveStress, NEMLObject, std::shared_ptr<EffectiveStress>>(m, "EffectiveStress")
      .def("effective",
           [](EffectiveStress & m, py::array_t<double, py::array::c_style> s) -> double
           {
            double res;
            int ier = m.effective(arr2ptr<double>(s), res);
            py_error(ier);
            return res;
           }, "The effective stress measure.")
      .def("deffective",
           [](EffectiveStress & m, py::array_t<double, py::array::c_style> s) -> py::array_t<double>
           {
            auto res = alloc_vec<double>(6);
            int ier = m.deffective(arr2ptr<double>(s), arr2ptr<double>(res));
            py_error(ier);
            return res;
           }, "The derivative of the effective stress measure wrt stress.")
      ;

  py::class_<VonMisesEffectiveStress, EffectiveStress, std::shared_ptr<VonMisesEffectiveStress>>(m, "VonMisesEffectiveStress")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<VonMisesEffectiveStress>(args, kwargs, {});
        }))
      ;

  py::class_<HuddlestonEffectiveStress, EffectiveStress, std::shared_ptr<HuddlestonEffectiveStress>>(m, "HuddlestonEffectiveStress")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<HuddlestonEffectiveStress>(args, kwargs, {"b"});
        }))
      ;

  py::class_<MaxPrincipalEffectiveStress, EffectiveStress, std::shared_ptr<MaxPrincipalEffectiveStress>>(m, "MaxPrincipalEffectiveStress")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<MaxPrincipalEffectiveStress>(args, kwargs, {});
        }))
      ;

  py::class_<MaxSeveralEffectiveStress, EffectiveStress, std::shared_ptr<MaxSeveralEffectiveStress>>(m, "MaxSeveralEffectiveStress")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<MaxSeveralEffectiveStress>(args, kwargs, {"measures"});
        }))
      ;

  py::class_<SumSeveralEffectiveStress, EffectiveStress, std::shared_ptr<SumSeveralEffectiveStress>>(m, "SumSeveralEffectiveStress")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<SumSeveralEffectiveStress>(args, kwargs, {"measures", "weights"});
        }))
      ;

  py::class_<ModularCreepDamageModel_sd, NEMLScalarDamagedModel_sd, std::shared_ptr<ModularCreepDamageModel_sd>>(m, "ModularCreepDamageModel_sd")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ModularCreepDamageModel_sd>(args, kwargs,
                                                                    {"elastic",
                                                                    "A", "xi",
                                                                    "phi", "estress", 
                                                                    "base"});
        }))
      ;

  py::class_<LarsonMillerCreepDamageModel_sd, NEMLScalarDamagedModel_sd, std::shared_ptr<LarsonMillerCreepDamageModel_sd>>(m, "LarsonMillerCreepDamageModel_sd")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<LarsonMillerCreepDamageModel_sd>(args, kwargs,
                                                                    {"elastic",
                                                                    "lmr",
                                                                    "estress", 
                                                                    "base"});
        }))
      ;

  py::class_<NEMLStandardScalarDamagedModel_sd, NEMLScalarDamagedModel_sd, std::shared_ptr<NEMLStandardScalarDamagedModel_sd>>(m, "NEMLStandardScalarDamagedModel_sd")
      .def("f",
           [](NEMLStandardScalarDamagedModel_sd & m, py::array_t<double, py::array::c_style> s_np1, double d_np1, double T_np1) -> double
           {
            double fv;
            int ier = m.f(arr2ptr<double>(s_np1), d_np1, T_np1, fv);
            py_error(ier);
            
            return fv;
           }, "The damage evolution function.")
      .def("df_ds",
           [](NEMLStandardScalarDamagedModel_sd & m, py::array_t<double, py::array::c_style> s_np1, double d_np1, double T_np1) -> py::array_t<double>
           {
            auto dfv = alloc_vec<double>(6);
            int ier = m.df_ds(arr2ptr<double>(s_np1), d_np1, T_np1, arr2ptr<double>(dfv));
            py_error(ier);

            return dfv;
           }, "The derivative of the damage function wrt. stress.")
      .def("df_dd",
           [](NEMLStandardScalarDamagedModel_sd & m, py::array_t<double, py::array::c_style> s_np1, double d_np1, double T_np1) -> double
           {
            double dfv;
            int ier = m.df_dd(arr2ptr<double>(s_np1), d_np1, T_np1, dfv);
            py_error(ier);
            
            return dfv;
           }, "The derivative of the damage function wrt. damage")

      ;

  py::class_<NEMLPowerLawDamagedModel_sd, NEMLStandardScalarDamagedModel_sd, std::shared_ptr<NEMLPowerLawDamagedModel_sd>>(m, "NEMLPowerLawDamagedModel_sd")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<NEMLPowerLawDamagedModel_sd>(args, kwargs, 
                                                                   {"elastic",
                                                                   "A", "a",
                                                                   "base"});
        }))
      ;

  py::class_<NEMLExponentialWorkDamagedModel_sd, NEMLStandardScalarDamagedModel_sd, std::shared_ptr<NEMLExponentialWorkDamagedModel_sd>>(m, "NEMLExponentialWorkDamagedModel_sd")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<NEMLExponentialWorkDamagedModel_sd>(
              args, kwargs, {"elastic", "W0", "k0", "af", "base"});
        }))
      ;
}

} //  namespace neml
