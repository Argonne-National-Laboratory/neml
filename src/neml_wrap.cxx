#include "neml.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_MODULE(neml, m) {
  m.doc() = "Base class for all material models.";
  
  py::class_<NEMLModel, std::shared_ptr<NEMLModel>>(m, "NEMLModel")
      .def_property_readonly("nstore", &NEMLModel::nstore, "Number of variables the program needs to store.")
      .def("init_store",
           [](NEMLModel & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nstore());
            int ier = m.init_store(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize stored variables.")

      .def_property_readonly("nhist", &NEMLModel::nhist, "Number of actual history variables.")
      .def("init_hist",
           [](NEMLModel & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nhist());
            int ier = m.init_hist(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize history variables.")

      .def("update_ldF",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> F_np1, py::array_t<double, py::array::c_style> F_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n, double u_n, double p_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>, double, double>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);
            double u_np1, p_np1;

            int ier = m.update_ldF(arr2ptr<double>(F_np1), arr2ptr<double>(F_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1), u_np1, u_n, p_np1, p_n);
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1, u_np1, p_np1);

           }, "Large deformation update through the deformation gradient.")

      .def("update_ldI",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> l_inc, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n, double u_n, double p_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>, double, double>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);
            double u_np1, p_np1;

            int ier = m.update_ldI(arr2ptr<double>(l_inc), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1), u_np1, u_n, p_np1, p_n);
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1, p_np1, p_n);

           }, "Large deformation incremental update through the spatial velocity gradient.")

      .def("update_sd",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n, double u_n, double p_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>, double, double>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);
            double u_np1, p_np1;

            int ier = m.update_sd(arr2ptr<double>(e_np1), arr2ptr<double>(e_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1), u_np1, u_n, p_np1, p_n);
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1, u_np1, p_np1);

           }, "Small deformation update.")

      .def("alpha", &NEMLModel::alpha)
      .def("elastic_strains",
           [](NEMLModel_sd & m, py::array_t<double, py::array::c_style> s_np1, double T_np1, py::array_t<double, py::array::c_style> h_np1) -> py::array_t<double>
           {
            auto e_np1 = alloc_vec<double>(6);

            int ier = m.elastic_strains(
                arr2ptr<double>(s_np1), T_np1,
                arr2ptr<double>(h_np1), arr2ptr<double>(e_np1));
            py_error(ier);

            return e_np1;

           }, "Calculate the elastic strains.")
      ;

  py::class_<NEMLModel_ldF, std::shared_ptr<NEMLModel_ldF>>(m, "NEMLModel_ldF", py::base<NEMLModel>())
      ;

  py::class_<NEMLModel_ldI, std::shared_ptr<NEMLModel_ldI>>(m, "NEMLModel_ldI", py::base<NEMLModel>())
      ;

  py::class_<NEMLModel_sd, std::shared_ptr<NEMLModel_sd>>(m, "NEMLModel_sd", py::base<NEMLModel>())
      ;

  py::class_<SmallStrainElasticity, std::shared_ptr<SmallStrainElasticity>>(m, "SmallStrainElasticity", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<LinearElasticModel>, double>(), 
           py::arg("elastic"), py::arg("alpha") = 0.0)
      .def(py::init<std::shared_ptr<LinearElasticModel>, 
           std::shared_ptr<Interpolate>>(),
           py::arg("elastic"), py::arg("alpha") = nullptr)
      .def_property_readonly("elastic", &SmallStrainElasticity::elastic)
      ;

  py::class_<SSPPTrialState>(m, "SSPPTrialState", py::base<TrialState>())
      ;

  py::class_<SSRIPTrialState>(m, "SSRIPTrialState", py::base<TrialState>())
      ;

  py::class_<SSCPTrialState>(m, "SSCPTrialState", py::base<TrialState>())
      ;

  py::class_<GITrialState>(m, "GITrialState", py::base<TrialState>())
      ;

  py::class_<SmallStrainPerfectPlasticity, std::shared_ptr<SmallStrainPerfectPlasticity>>(m, "SmallStrainPerfectPlasticity", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<LinearElasticModel>, 
           std::shared_ptr<YieldSurface>, 
           double, double,
           double, int , bool, int>(),
           py::arg("elastic"), py::arg("flow"), py::arg("ys"), 
           py::arg("alpha") = 0.0,
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false, py::arg("max_divide") = 8)
 
      .def(py::init<std::shared_ptr<LinearElasticModel>, 
           std::shared_ptr<YieldSurface>, 
           std::shared_ptr<Interpolate>,
           std::shared_ptr<Interpolate>,
           double, int , bool, int>(),
           py::arg("elastic"), py::arg("flow"), py::arg("ys"), 
           py::arg("alpha") = nullptr,
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false, py::arg("max_divide") = 8)

      .def("ys", &SmallStrainPerfectPlasticity::ys)
      .def_property_readonly("elastic", &SmallStrainPerfectPlasticity::elastic)

      .def("make_trial_state",
           [](SmallStrainPerfectPlasticity & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::unique_ptr<SSPPTrialState>
           {
              std::unique_ptr<SSPPTrialState> ts(new SSPPTrialState);
              int ier = m.make_trial_state(arr2ptr<double>(e_np1),
                                          arr2ptr<double>(e_n),
                                          T_np1, T_n,
                                          t_np1, t_n,
                                          arr2ptr<double>(s_n),
                                          arr2ptr<double>(h_n),
                                          *ts);
              py_error(ier);

              return ts;
           }, "Setup trial state for solve.")

      // Remove if/when pybind11 supports multiple inheritance
      .def_property_readonly("nparams", &SmallStrainPerfectPlasticity::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](SmallStrainPerfectPlasticity & m, SSPPTrialState & ts) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x), &ts);
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](SmallStrainPerfectPlasticity & m, py::array_t<double, py::array::c_style> x, SSPPTrialState & ts) -> std::tuple<py::array_t<double>, py::array_t<double>>
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

  py::class_<SmallStrainRateIndependentPlasticity, std::shared_ptr<SmallStrainRateIndependentPlasticity>>(m, "SmallStrainRateIndependentPlasticity", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<LinearElasticModel>, std::shared_ptr<RateIndependentFlowRule>, double, double, int , bool, double, bool>(),
           py::arg("elastic"), py::arg("flow"),
           py::arg("alpha") = 0.0,
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false, py::arg("kttol") = 1.0e-2,
           py::arg("check_kt") = false)

      .def(py::init<std::shared_ptr<LinearElasticModel>, std::shared_ptr<RateIndependentFlowRule>, std::shared_ptr<Interpolate>, double, int , bool, double, bool>(),
           py::arg("elastic"), py::arg("flow"),
           py::arg("alpha") = nullptr,
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false, py::arg("kttol") = 1.0e-2,
           py::arg("check_kt") = false)

      .def_property_readonly("elastic", &SmallStrainRateIndependentPlasticity::elastic)

      .def("make_trial_state",
           [](SmallStrainRateIndependentPlasticity & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::unique_ptr<SSRIPTrialState>
           {
              std::unique_ptr<SSRIPTrialState> ts(new SSRIPTrialState);
              int ier = m.make_trial_state(arr2ptr<double>(e_np1),
                                          arr2ptr<double>(e_n),
                                          T_np1, T_n,
                                          t_np1, t_n,
                                          arr2ptr<double>(s_n),
                                          arr2ptr<double>(h_n),
                                          *ts);
              py_error(ier);

              return ts;
           }, "Setup trial state for solve.")

      // Remove if/when pybind11 supports multiple inheritance
      .def_property_readonly("nparams", &SmallStrainRateIndependentPlasticity::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](SmallStrainRateIndependentPlasticity & m, SSRIPTrialState & ts) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x), &ts);
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](SmallStrainRateIndependentPlasticity & m, py::array_t<double, py::array::c_style> x, SSRIPTrialState & ts) -> std::tuple<py::array_t<double>, py::array_t<double>>
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

  py::class_<SmallStrainCreepPlasticity, std::shared_ptr<SmallStrainCreepPlasticity>>(m, "SmallStrainCreepPlasticity", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<NEMLModel_sd>, std::shared_ptr<CreepModel>, double, double, int , bool, double>(),
           py::arg("plastic"), py::arg("creep"),
           py::arg("alpha") = 0.0,
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false, py::arg("sf") = 1.0e6)

      .def(py::init<std::shared_ptr<NEMLModel_sd>, std::shared_ptr<CreepModel>, std::shared_ptr<Interpolate>, double, int , bool, double>(),
           py::arg("plastic"), py::arg("creep"),
           py::arg("alpha") = nullptr,
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false, py::arg("sf") = 1.0e6)

      .def("make_trial_state",
           [](SmallStrainCreepPlasticity & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::unique_ptr<SSCPTrialState>
           {
              std::unique_ptr<SSCPTrialState> ts(new SSCPTrialState);
              int ier = m.make_trial_state(arr2ptr<double>(e_np1),
                                          arr2ptr<double>(e_n),
                                          T_np1, T_n,
                                          t_np1, t_n,
                                          arr2ptr<double>(s_n),
                                          arr2ptr<double>(h_n),
                                          *ts);
              py_error(ier);

              return ts;
           }, "Setup trial state for solve.")

      // Remove if/when pybind11 supports multiple inheritance
      .def_property_readonly("nparams", &SmallStrainCreepPlasticity::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](SmallStrainCreepPlasticity & m, SSRIPTrialState & ts) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x), &ts);
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](SmallStrainCreepPlasticity & m, py::array_t<double, py::array::c_style> x, SSCPTrialState & ts) -> std::tuple<py::array_t<double>, py::array_t<double>>
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



  py::class_<GeneralIntegrator, std::shared_ptr<GeneralIntegrator>>(m, "GeneralIntegrator", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<GeneralFlowRule>, double, double, int , bool, int>(),
           py::arg("rule"),
           py::arg("alpha") = 0.0,
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false,
           py::arg("max_divide") = 6)

      .def(py::init<std::shared_ptr<GeneralFlowRule>, std::shared_ptr<Interpolate>, double, int , bool, int>(),
           py::arg("rule"),
           py::arg("alpha") = nullptr,
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false,
           py::arg("max_divide") = 6)
  
      .def("make_trial_state",
           [](GeneralIntegrator & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::unique_ptr<GITrialState>
           {
              std::unique_ptr<GITrialState> ts(new GITrialState);
              int ier = m.make_trial_state(arr2ptr<double>(e_np1),
                                          arr2ptr<double>(e_n),
                                          T_np1, T_n,
                                          t_np1, t_n,
                                          arr2ptr<double>(s_n),
                                          arr2ptr<double>(h_n),
                                          *ts);
              py_error(ier);
              return ts;
           }, "Setup trial state for solve.")

      // Remove if/when pybind11 supports multiple inheritance
      .def_property_readonly("nparams", &GeneralIntegrator::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](GeneralIntegrator & m, GITrialState & ts) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x), &ts);
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](GeneralIntegrator & m, py::array_t<double, py::array::c_style> x, GITrialState & ts) -> std::tuple<py::array_t<double>, py::array_t<double>>
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

  py::class_<KMRegimeModel, std::shared_ptr<KMRegimeModel>>(m, "KMRegimeModel", py::base<NEMLModel_sd>())
      .def(py::init<std::vector<std::shared_ptr<NEMLModel_sd>>,
           std::vector<double>, std::shared_ptr<LinearElasticModel>,
           double, double, double, double>(),
           py::arg("models"), py::arg("gs"), py::arg("emodel"),
           py::arg("kboltz"), py::arg("b"), py::arg("eps0"),
           py::arg("alpha") = 0.0)
      .def(py::init<std::vector<std::shared_ptr<NEMLModel_sd>>,
           std::vector<double>, std::shared_ptr<LinearElasticModel>,
           double, double, double, std::shared_ptr<Interpolate>>(),
           py::arg("models"), py::arg("gs"), py::arg("emodel"),
           py::arg("kboltz"), py::arg("b"), py::arg("eps0"),
           py::arg("alpha") = nullptr)
      ;
}

} // namespace neml
