#include "volumetric.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(volumetric) {
  py::module m("volumetric", "Volumetric small strain material models.");

  py::class_<VolumetricModel, std::shared_ptr<VolumetricModel>>(m, "VolumetricModel")
      .def("update",
           [](const VolumetricModel & m, py::array_t<double, py::array::c_style> e_np1, py::array_t<double, py::array::c_style> e_n, double T_np1, double T_n, double t_np1, double t_n, py::array_t<double, py::array::c_style> s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
           {
            auto h_np1 = alloc_vec<double>(m.nhist());
            auto s_np1 = alloc_vec<double>(6);
            auto A_np1 = alloc_mat<double>(6, 6);
            
            int ier = m.update(arr2ptr<double>(e_np1), arr2ptr<double>(e_n), T_np1, T_n, t_np1, t_n, arr2ptr<double>(s_np1), arr2ptr<double>(s_n), arr2ptr<double>(h_np1), arr2ptr<double>(h_n), arr2ptr<double>(A_np1));
            raise_py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1);
           }, "Vector volumetric update")

      .def_property_readonly("nhist", &VolumetricModel::nhist)

      .def("init_hist",
           [](const VolumetricModel & m) -> py::array_t<double>
           {
            auto h = alloc_vec<double>(m.nhist());
            m.init_hist(arr2ptr<double>(h));
            return h;
           }, "Initialize history.")

      .def("update_mean", 
           [](const VolumetricModel & m, double e_np1, double e_n, double T_np1, double T_n, double t_np1, double t_n, double s_n, py::array_t<double, py::array::c_style> h_n) -> std::tuple<double, py::array_t<double>, double>
           {
            double s_np1, A_np1;
            auto h_np1 = alloc_vec<double>(m.nhist());

            int ier =  m.update_mean(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n, arr2ptr<double>(h_np1), arr2ptr<double>(h_n), A_np1);
            raise_py_error(ier);

            return std::make_tuple(s_np1, h_np1, A_np1);
           }, "Scalar volumetric strain/mean stress update")
      ;

  py::class_<VModelConstantK, std::shared_ptr<VModelConstantK>>(m, "VModelConstantK", py::base<VolumetricModel>())
      .def(py::init<double>())
      .def_property_readonly("K", &VModelConstantK::K)
      ;

  return m.ptr();
}

} // namespace neml
