#include "volumetric.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

namespace neml {

PYBIND11_PLUGIN(volumetric) {
  py::module m("volumetric", "Volumetric small strain material models.");

  py::class_<VolumetricModel>(m, "VolumetricModel")
      .def("update",
           [](const VolumetricModel & m, py::array_t<double, py::array::c_style> e_inc, double T_np1, double T_inc, double t_np1, double t_inc, py::array_t<double, py::array::c_style> h_n, py::array_t<double, py::array::c_style> s_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
           {
            auto h_np1 = alloc_vec<double>(m.nhist());
            auto s_np1 = alloc_vec<double>(6);
            auto A_np1 = alloc_mat<double>(6, 6);
            
            int ier = m.update(arr2ptr<double>(e_inc), T_np1, T_inc, t_np1, t_inc, arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(A_np1));
            raise_py_error(ier);

            return std::make_tuple(h_np1, s_np1, A_np1);
           }, "Vector volumetric update")
      .def_property_readonly("nhist", &VolumetricModel::nhist)
      .def("update_mean", 
           [](const VolumetricModel & m, double e_inc, double T_np1, double T_inc, double t_np1, double t_inc, py::array_t<double, py::array::c_style> h_n, double s_n) -> std::tuple<py::array_t<double>, double, double>
           {
            double s_np1, A_np1;
            auto h_np1 = alloc_vec<double>(m.nhist());

            int ier =  m.update_mean(e_inc, T_np1, T_inc, t_np1, t_inc, arr2ptr<double>(h_np1), arr2ptr<double>(h_n), s_np1, s_n, A_np1);
            raise_py_error(ier);

            return std::make_tuple(h_np1, s_np1, A_np1);
           }, "Scalar volumetric strain/mean stress update")
      ;

  py::class_<VModelConstantK>(m, "VModelConstantK", py::base<VolumetricModel>())
      .def(py::init<double>())
      .def_property_readonly("K", &VModelConstantK::K)
      ;

  return m.ptr();
}

} // namespace neml
