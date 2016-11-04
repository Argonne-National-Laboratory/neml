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
           [](const NEMLModel & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nstore());
            int ier = m.init_store(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize stored variables.")

      .def_property_readonly("nhist", &NEMLModel::nhist, "Number of actual history variables.")
      .def("init_store",
           [](const NEMLModel & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nhist());
            int ier = m.init_hist(arr2ptr<double>(h));
            py_error(ier);
            return h;
           }, "Initialize history variables.")

      .def("update_ldF",
           [](const NEMLModel & m, py::array_t<double, py::array::c_style> F_np1, py::array_t<double, py::array::c_style> F_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);

            int ier = m.update_ldF(arr2ptr<double>(F_np1), arr2ptr<double>(F_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1));
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1);

           }, "Large deformation update through the deformation gradient.")

      .def("update_ldI",
           [](const NEMLModel & m, py::array_t<double, py::array::c_style> l_inc, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
           {
            auto s_np1 = alloc_vec<double>(6);
            auto h_np1 = alloc_vec<double>(m.nstore());
            auto A_np1 = alloc_mat<double>(6,6);

            int ier = m.update_ldI(arr2ptr<double>(l_inc), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1));
            py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1);

           }, "Large deformation incremental update through the spatial velocity gradient.")

      .def("update_sd",
           [](const NEMLModel & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
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

  py::class_<SplitModel_sd, std::shared_ptr<SplitModel_sd>>(m, "SplitModel_sd", py::base<NEMLModel_sd>())
      .def(py::init<std::shared_ptr<VolumetricModel>, std::shared_ptr<DeviatoricModel>>())
      ;

  return m.ptr();
}


} // namespace neml
