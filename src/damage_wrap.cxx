#include "pyhelp.h" // include first to avoid annoying redef warning

#include "damage.h"

#include "nemlerror.h"
#include "parse.h"

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
      .def("init_damage", &NEMLDamagedModel_sd::init_damage)
      .def("populate_damage", &NEMLDamagedModel_sd::populate_damage)
      ;

  py::class_<SDTrialState, TrialState>(m, "SDTrialState")
      ;

  py::class_<NEMLScalarDamagedModel_sd, NEMLDamagedModel_sd, Solvable, std::shared_ptr<NEMLScalarDamagedModel_sd>>(m, "NEMLScalarDamagedModel_sd")
      PICKLEABLE(NEMLScalarDamagedModel_sd)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<NEMLScalarDamagedModel_sd>(args, kwargs,
                                                              {"elastic", "base",
                                                              "damage"});
        }))
      .def("make_trial_state",
           [](NEMLScalarDamagedModel_sd & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n, double u_n, double p_n) -> std::unique_ptr<SDTrialState>
           {
            std::unique_ptr<SDTrialState> ts(new SDTrialState);
            m.make_trial_state(arr2ptr<double>(e_np1),
                                         arr2ptr<double>(e_n),
                                         T_np1, T_n, t_np1, t_n,
                                         arr2ptr<double>(s_n),
                                         arr2ptr<double>(h_n),
                                         u_n, p_n, *ts);

            return ts;
           }, "Make a trial state, mostly for testing.")
      ;

  py::class_<ScalarDamage, NEMLObject, std::shared_ptr<ScalarDamage>>(m,
                                                                      "ScalarDamage")
      .def("damage",
           [](ScalarDamage & m, double d_np1, double d_n, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> s_n, double T_np1, double T_n, double t_np1, double t_n) -> double
           {
            double damage;
            m.damage(d_np1, d_n, arr2ptr<double>(e_np1), arr2ptr<double>(e_n),
                     arr2ptr<double>(s_np1), arr2ptr<double>(s_n),
                     T_np1, T_n,
                     t_np1, t_n, &damage);
            return damage;
           }, "The damage evolution equation.")
      .def("ddamage_dd",
           [](ScalarDamage & m, double d_np1, double d_n, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> s_n, double T_np1, double T_n, double t_np1, double t_n) -> double
           {
            double ddamage;
            m.ddamage_dd(d_np1, d_n, arr2ptr<double>(e_np1), arr2ptr<double>(e_n),
                     arr2ptr<double>(s_np1), arr2ptr<double>(s_n),
                     T_np1, T_n,
                     t_np1, t_n, &ddamage);
            return ddamage;
           }, "The derivative of the damage evolution equation wrt. damage.")
      .def("ddamage_de",
           [](ScalarDamage & m, double d_np1, double d_n, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> s_n, double T_np1, double T_n, double t_np1, double t_n) -> py::array_t<double>
           {
            auto ddamage = alloc_vec<double>(6);
            m.ddamage_de(d_np1, d_n, arr2ptr<double>(e_np1), arr2ptr<double>(e_n),
                     arr2ptr<double>(s_np1), arr2ptr<double>(s_n),
                     T_np1, T_n,
                     t_np1, t_n, arr2ptr<double>(ddamage));
            return ddamage;
           }, "The derivative of the damage evolution equation wrt. strain.") 
      .def("ddamage_ds",
           [](ScalarDamage & m, double d_np1, double d_n, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_np1, py::array_t<double, py::array::c_style> s_n, double T_np1, double T_n, double t_np1, double t_n) -> py::array_t<double>
           {
            auto ddamage = alloc_vec<double>(6);
            m.ddamage_ds(d_np1, d_n, arr2ptr<double>(e_np1), arr2ptr<double>(e_n),
                     arr2ptr<double>(s_np1), arr2ptr<double>(s_n),
                     T_np1, T_n,
                     t_np1, t_n, arr2ptr<double>(ddamage));
            return ddamage;
           }, "The derivative of the damage evolution equation wrt. stress.")
      .def_property_readonly("d_init", &ScalarDamage::d_init)
      ;

  py::class_<ScalarDamageRate, ScalarDamage, std::shared_ptr<ScalarDamageRate>>(m, "ScalarDamageRate")
      .def("damage_rate",
           [](ScalarDamageRate & m, double d, py::array_t<double, py::array::c_style> e, py::array_t<double, py::array::c_style> s, double T, double t) -> double
           {
            double drate;
            m.damage_rate(d, arr2ptr<double>(e), arr2ptr<double>(s), T, t, &drate);
            return drate;
           }, "The damage rate")
      .def("ddamage_rate_dd",
           [](ScalarDamageRate & m, double d, py::array_t<double, py::array::c_style> e, py::array_t<double, py::array::c_style> s, double T, double t) -> double
           {
            double drate;
            m.ddamage_rate_dd(d, arr2ptr<double>(e), arr2ptr<double>(s), T, t, &drate);
            return drate;
           }, "The derivative of the damage rate with respect to the scalar damage") 
      .def("ddamage_rate_de",
           [](ScalarDamageRate & m, double d, py::array_t<double, py::array::c_style> e, py::array_t<double, py::array::c_style> s, double T, double t) -> py::array_t<double>
           {
            auto ddamage = alloc_vec<double>(6);
            m.ddamage_rate_de(d, arr2ptr<double>(e), arr2ptr<double>(s), T, t, arr2ptr<double>(ddamage));
            return ddamage;
           }, "The derivative of the damage rate with respect to the strain")
      .def("ddamage_rate_ds",
           [](ScalarDamageRate & m, double d, py::array_t<double, py::array::c_style> e, py::array_t<double, py::array::c_style> s, double T, double t) -> py::array_t<double>
           {
            auto ddamage = alloc_vec<double>(6);
            m.ddamage_rate_ds(d, arr2ptr<double>(e), arr2ptr<double>(s), T, t, arr2ptr<double>(ddamage));
            return ddamage;
           }, "The derivative of the damage rate with respect to the stress") 
      ;

  py::class_<CombinedDamage, ScalarDamage, std::shared_ptr<CombinedDamage>>(m, "CombinedDamage")
      PICKLEABLE(CombinedDamage)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<CombinedDamage>(args, kwargs,
                                                              {"elastic",
                                                              "models"});
        }))
      ;

  py::class_<ClassicalCreepDamage, ScalarDamageRate, std::shared_ptr<ClassicalCreepDamage>>(m, "ClassicalCreepDamage")
      PICKLEABLE(ClassicalCreepDamage)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ClassicalCreepDamage>(args, kwargs,
                                                                    {"elastic",
                                                                    "A", "xi",
                                                                    "phi"});
        }))
      ;

  py::class_<WorkDamage, ScalarDamage,
      std::shared_ptr<WorkDamage>>(m, "WorkDamage")
      PICKLEABLE(WorkDamage)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<WorkDamage>(args, kwargs,
                                                                    {"elastic",
                                                                    "Wcrit", "n"});
        }))
      ;

  py::class_<EffectiveStress, NEMLObject, std::shared_ptr<EffectiveStress>>(m, "EffectiveStress")
      .def("effective",
           [](EffectiveStress & m, py::array_t<double, py::array::c_style> s) -> double
           {
            double res;
            m.effective(arr2ptr<double>(s), res);
            return res;
           }, "The effective stress measure.")
      .def("deffective",
           [](EffectiveStress & m, py::array_t<double, py::array::c_style> s) -> py::array_t<double>
           {
            auto res = alloc_vec<double>(6);
            m.deffective(arr2ptr<double>(s), arr2ptr<double>(res));
            return res;
           }, "The derivative of the effective stress measure wrt stress.")
      ;

  py::class_<VonMisesEffectiveStress, EffectiveStress, std::shared_ptr<VonMisesEffectiveStress>>(m, "VonMisesEffectiveStress")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<VonMisesEffectiveStress>(args, kwargs, {});
        }))
      ;

  py::class_<MeanEffectiveStress, EffectiveStress, std::shared_ptr<MeanEffectiveStress>>(m, "MeanEffectiveStress")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<MeanEffectiveStress>(args, kwargs, {});
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

  py::class_<ModularCreepDamage, ScalarDamageRate, std::shared_ptr<ModularCreepDamage>>(m, "ModularCreepDamage")
      PICKLEABLE(ModularCreepDamage)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ModularCreepDamage>(args, kwargs,
                                                                    {"elastic",
                                                                    "A", "xi",
                                                                    "phi", "estress"});
        }))
      ;

  py::class_<LarsonMillerCreepDamage, ScalarDamageRate, std::shared_ptr<LarsonMillerCreepDamage>>(m, "LarsonMillerCreepDamage")
      PICKLEABLE(LarsonMillerCreepDamage)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<LarsonMillerCreepDamage>(args, kwargs,
                                                                    {"elastic",
                                                                    "lmr",
                                                                    "estress"});
        }))
      ;

  py::class_<StandardScalarDamage, ScalarDamage, std::shared_ptr<StandardScalarDamage>>(m, "StandardScalarDamage")
      .def("f",
           [](StandardScalarDamage & m, py::array_t<double, py::array::c_style> s_np1, double d_np1, double T_np1) -> double
           {
            double fv;
            m.f(arr2ptr<double>(s_np1), d_np1, T_np1, fv);
            
            return fv;
           }, "The damage evolution function.")
      .def("df_ds",
           [](StandardScalarDamage & m, py::array_t<double, py::array::c_style> s_np1, double d_np1, double T_np1) -> py::array_t<double>
           {
            auto dfv = alloc_vec<double>(6);
            m.df_ds(arr2ptr<double>(s_np1), d_np1, T_np1, arr2ptr<double>(dfv));

            return dfv;
           }, "The derivative of the damage function wrt. stress.")
      .def("df_dd",
           [](StandardScalarDamage & m, py::array_t<double, py::array::c_style> s_np1, double d_np1, double T_np1) -> double
           {
            double dfv;
            m.df_dd(arr2ptr<double>(s_np1), d_np1, T_np1, dfv);
            
            return dfv;
           }, "The derivative of the damage function wrt. damage")

      ;

  py::class_<PowerLawDamage, StandardScalarDamage, std::shared_ptr<PowerLawDamage>>(m, "PowerLawDamage")
      PICKLEABLE(PowerLawDamage)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<PowerLawDamage>(args, kwargs, 
                                                                   {"elastic",
                                                                   "A", "a"});
        }))
      ;

  py::class_<ExponentialWorkDamage, StandardScalarDamage, std::shared_ptr<ExponentialWorkDamage>>(m, "ExponentialWorkDamage")
      PICKLEABLE(ExponentialWorkDamage)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ExponentialWorkDamage>(
              args, kwargs, {"elastic", "W0", "k0", "af"});
        }))
      ;
}

} //  namespace neml
