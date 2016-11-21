#include "neml.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(neml) {
  py::module m("neml", "Base class material models.");
  
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
      .def("init_store",
           [](NEMLModel & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nhist());
            int ier = m.init_hist(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize history variables.")

      .def("update_ldF",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> F_np1, py::array_t<double, py::array::c_style> F_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);

            int ier = m.update_ldF(arr2ptr<double>(F_np1), arr2ptr<double>(F_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1));
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1);

           }, "Large deformation update through the deformation gradient.")

      .def("update_ldI",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> l_inc, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);

            int ier = m.update_ldI(arr2ptr<double>(l_inc), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1));
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1);

           }, "Large deformation incremental update through the spatial velocity gradient.")

      .def("update_sd",
           [](NEMLModel & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);

            int ier = m.update_sd(arr2ptr<double>(e_np1), arr2ptr<double>(e_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1));
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1);

           }, "Small deformation update.")
      ;

  py::class_<NEMLModel_ldF, std::shared_ptr<NEMLModel_ldF>>(m, "NEMLModel_ldF", py::base<NEMLModel>())
      ;

  py::class_<NEMLModel_ldI, std::shared_ptr<NEMLModel_ldI>>(m, "NEMLModel_ldI", py::base<NEMLModel>())
      ;

  py::class_<NEMLModel_sd, std::shared_ptr<NEMLModel_sd>>(m, "NEMLModel_sd", py::base<NEMLModel>())
      ;

  py::class_<SmallStrainElasticity, std::shared_ptr<SmallStrainElasticity>>(m, "SmallStrainElasticity", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<LinearElasticModel>>(), 
           py::arg("elastic"))
      ;

  py::class_<SmallStrainRateIndependentPlasticity, std::shared_ptr<SmallStrainRateIndependentPlasticity>>(m, "SmallStrainRateIndependentPlasticity", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<LinearElasticModel>, std::shared_ptr<RateIndependentFlowRule>, double, int , bool, double, bool>(),
           py::arg("elastic"), py::arg("flow"), 
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false, py::arg("kttol") = 1.0e-2,
           py::arg("check_kt") = true)
  
      .def("set_trial_state",
           [](SmallStrainRateIndependentPlasticity & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> h_n, double T_np1, double t_np1, double t_n) -> void
           {
              int ier = m.set_trial_state(arr2ptr<double>(e_np1), arr2ptr<double>(h_n), T_np1, t_np1, t_n);
              py_error(ier);
           }, "Setup trial state for solve.")

      // Remove if/when pybind11 supports multiple inheritance
      .def_property_readonly("nparams", &SmallStrainRateIndependentPlasticity::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](SmallStrainRateIndependentPlasticity & m) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x));
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](SmallStrainRateIndependentPlasticity & m, py::array_t<double, py::array::c_style> x) -> std::tuple<py::array_t<double>, py::array_t<double>>
           {
            auto R = alloc_vec<double>(m.nparams());
            auto J = alloc_mat<double>(m.nparams(), m.nparams());
            
            int ier = m.RJ(arr2ptr<double>(x), arr2ptr<double>(R), arr2ptr<double>(J));
            py_error(ier);

            return std::make_tuple(R, J);
           }, "Residual and jacobian.")
      // End remove block
      ;

  py::class_<GeneralIntegrator, std::shared_ptr<GeneralIntegrator>>(m, "GeneralIntegrator", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<GeneralFlowRule>, double, int , bool>(),
           py::arg("rule"),
           py::arg("tol") = 1.0e-8, py::arg("miter") = 50, 
           py::arg("verbose") = false)
  
      .def("set_trial_state",
           [](GeneralIntegrator & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n, double T_np1, double T_n, double t_np1, double t_n) -> void
           {
              int ier = m.set_trial_state(arr2ptr<double>(e_np1), arr2ptr<double>(e_n), arr2ptr<double>(s_n), arr2ptr<double>(h_n), T_np1, T_n, t_np1, t_n);
              py_error(ier);
           }, "Setup trial state for solve.")

      // Remove if/when pybind11 supports multiple inheritance
      .def_property_readonly("nparams", &GeneralIntegrator::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](GeneralIntegrator & m) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x));
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](GeneralIntegrator & m, py::array_t<double, py::array::c_style> x) -> std::tuple<py::array_t<double>, py::array_t<double>>
           {
            auto R = alloc_vec<double>(m.nparams());
            auto J = alloc_mat<double>(m.nparams(), m.nparams());
            
            int ier = m.RJ(arr2ptr<double>(x), arr2ptr<double>(R), arr2ptr<double>(J));
            py_error(ier);

            return std::make_tuple(R, J);
           }, "Residual and jacobian.")
      // End remove block
      ;

  return m.ptr();
}


} // namespace neml
