#include "pyhelp.h" // include first to avoid annoying redef warning

#include "models.h"

#include "nemlerror.h"
#include "parse.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(models, m) {
  py::module::import("neml.objects");
  py::module::import("neml.solvers");

  m.doc() = "Base class for all material models.";
  
  py::class_<NEMLModel, HistoryNEMLObject, std::shared_ptr<NEMLModel>>(m, "NEMLModel")
      .def("save", &NEMLModel::save)
      .def("report_internal_variable_names", &NEMLModel::report_internal_variable_names)
      .def("update_sd",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n, double u_n, double p_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>, double, double>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);
            double u_np1, p_np1;

            m.update_sd(arr2ptr<double>(e_np1), arr2ptr<double>(e_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1), u_np1, u_n, p_np1, p_n);

            return std::make_tuple(s_np1, h_np1, A_np1, u_np1, p_np1);

           }, "Small deformation update.")
      .def("update_ld_inc",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> d_np1, py::array_t<double, py::array::c_style> d_n, py::array_t<double, py::array::c_style> w_np1, py::array_t<double, py::array::c_style> w_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n, double u_n, double p_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>, py::array_t<double>, double, double>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);
            auto B_np1 = alloc_mat<double>(6,3);
            double u_np1, p_np1;

            m.update_ld_inc(arr2ptr<double>(d_np1), arr2ptr<double>(d_n), arr2ptr<double>(w_np1), arr2ptr<double>(w_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1), arr2ptr<double>(B_np1), u_np1, u_n, p_np1, p_n);

            return std::make_tuple(s_np1, h_np1, A_np1, B_np1, u_np1, p_np1);

           }, "Large deformation incremental update.")

      .def("alpha", &NEMLModel::alpha)
      .def("elastic_strains",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> s_np1, double T_np1, py::array_t<double, py::array::c_style> h_np1) -> py::array_t<double>
           {
            auto e_np1 = alloc_vec<double>(6);

            m.elastic_strains(
                arr2ptr<double>(s_np1), T_np1,
                arr2ptr<double>(h_np1), arr2ptr<double>(e_np1));

            return e_np1;

           }, "Calculate the elastic strains.")
      .def("get_damage",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> h_np1) -> double
           {
            return m.get_damage(arr2ptr<double>(h_np1));
           }, "Provide the current value of damage, for models with damage")
      .def("should_del_element",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> h_np1) -> bool
           {
            return m.should_del_element(arr2ptr<double>(h_np1));
           }, "Report if element death should occur")
      .def("is_damage_model", &NEMLModel::is_damage_model, "Report if this model is a damage model")
      ;

  py::class_<NEMLModel_sd, NEMLModel, std::shared_ptr<NEMLModel_sd>>(m, "NEMLModel_sd")
      .def_property_readonly("elastic", &NEMLModel_sd::elastic)
      .def("set_elastic_model", &NEMLModel_sd::set_elastic_model)
      ;

  py::class_<NEMLModel_ldi, NEMLModel, std::shared_ptr<NEMLModel_ldi>>(m, "NEMLModel_ldi")
      ;

  py::class_<SubstepModel_sd, NEMLModel_sd, std::shared_ptr<SubstepModel_sd>>(m,
                                                                              "SubstepModel_sd")
      ;

  py::class_<SmallStrainElasticity, NEMLModel_sd, std::shared_ptr<SmallStrainElasticity>>(m, "SmallStrainElasticity")
      PICKLEABLE(SmallStrainElasticity)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<SmallStrainElasticity>(args, kwargs,
                                                             {"elastic"});
        }))
      ;

  py::class_<SSPPTrialState, TrialState>(m, "SSPPTrialState")
      ;

  py::class_<SSRIPTrialState, TrialState>(m, "SSRIPTrialState")
      ;

  py::class_<SSCPTrialState, TrialState>(m, "SSCPTrialState")
      ;

  py::class_<GITrialState, TrialState>(m, "GITrialState")
      ;

  py::class_<SmallStrainPerfectPlasticity, SubstepModel_sd, Solvable, std::shared_ptr<SmallStrainPerfectPlasticity>>(m, "SmallStrainPerfectPlasticity")
      PICKLEABLE(SmallStrainPerfectPlasticity)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<SmallStrainPerfectPlasticity>(args, kwargs,
                                                             {"elastic", "surface",
                                                             "ys"});
        }))

      .def("ys", &SmallStrainPerfectPlasticity::ys)

      .def("make_trial_state",
           [](SmallStrainPerfectPlasticity & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::unique_ptr<SSPPTrialState>
           {
              std::unique_ptr<SSPPTrialState> ts(new SSPPTrialState);
              m.make_trial_state(arr2ptr<double>(e_np1),
                                          arr2ptr<double>(e_n),
                                          T_np1, T_n,
                                          t_np1, t_n,
                                          arr2ptr<double>(s_n),
                                          arr2ptr<double>(h_n),
                                          *ts);
              return ts;
           }, "Setup trial state for solve.")
      ;

  py::class_<SmallStrainRateIndependentPlasticity, SubstepModel_sd, Solvable, std::shared_ptr<SmallStrainRateIndependentPlasticity>>(m, "SmallStrainRateIndependentPlasticity")
      PICKLEABLE(SmallStrainRateIndependentPlasticity)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<SmallStrainRateIndependentPlasticity>(args, kwargs,
                                                             {"elastic", "flow"});
        }))

      .def("make_trial_state",
           [](SmallStrainRateIndependentPlasticity & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::unique_ptr<SSRIPTrialState>
           {
              std::unique_ptr<SSRIPTrialState> ts(new SSRIPTrialState);
              m.make_trial_state(arr2ptr<double>(e_np1),
                                          arr2ptr<double>(e_n),
                                          T_np1, T_n,
                                          t_np1, t_n,
                                          arr2ptr<double>(s_n),
                                          arr2ptr<double>(h_n),
                                          *ts);
              return ts;
           }, "Setup trial state for solve.")
      ;

  py::class_<SmallStrainCreepPlasticity, NEMLModel_sd, Solvable, std::shared_ptr<SmallStrainCreepPlasticity>>(m, "SmallStrainCreepPlasticity")
      PICKLEABLE(SmallStrainCreepPlasticity)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<SmallStrainCreepPlasticity>(args, kwargs,
                                                                  {"elastic",
                                                                  "plastic",
                                                                  "creep"});
        }))
      .def("make_trial_state",
           [](SmallStrainCreepPlasticity & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::unique_ptr<SSCPTrialState>
           {
              std::unique_ptr<SSCPTrialState> ts(new SSCPTrialState);
              m.make_trial_state(arr2ptr<double>(e_np1),
                                          arr2ptr<double>(e_n),
                                          T_np1, T_n,
                                          t_np1, t_n,
                                          arr2ptr<double>(s_n),
                                          arr2ptr<double>(h_n),
                                          *ts);
              return ts;
           }, "Setup trial state for solve.")
      ;



  py::class_<GeneralIntegrator, SubstepModel_sd, Solvable, std::shared_ptr<GeneralIntegrator>>(m, "GeneralIntegrator")
      PICKLEABLE(GeneralIntegrator)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<GeneralIntegrator>(args, kwargs, 
                                                         {"elastic", "rule"});
        }))

      .def("make_trial_state",
           [](GeneralIntegrator & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::unique_ptr<GITrialState>
           {
              std::unique_ptr<GITrialState> ts(new GITrialState);
              m.make_trial_state(arr2ptr<double>(e_np1),
                                          arr2ptr<double>(e_n),
                                          T_np1, T_n,
                                          t_np1, t_n,
                                          arr2ptr<double>(s_n),
                                          arr2ptr<double>(h_n),
                                          *ts);
              return ts;
           }, "Setup trial state for solve.")
      ;

  py::class_<KMRegimeModel, NEMLModel_sd, std::shared_ptr<KMRegimeModel>>(m, "KMRegimeModel")
      PICKLEABLE(KMRegimeModel)
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<KMRegimeModel>(args, kwargs, 
                                                     {"elastic", "models",
                                                     "gs", "kboltz",
                                                     "b", "eps0"});
        }))
      ;
}

} // namespace neml
