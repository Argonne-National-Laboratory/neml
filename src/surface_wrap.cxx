#include "surface.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

namespace neml {

PYBIND11_PLUGIN(surface) {
  py::module m("surface", "Various yield surfaces and associated hardening rules.");

  py::class_<YieldSurface>(m, "YieldSurface")
      .def_property_readonly("nhist", &YieldSurface::nhist, "Number of history variables.")

      .def("f", 
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> double
           {
            double fv;

            int ier = m.f(arr2ptr<double>(s), arr2ptr<double>(h), T, fv);

            return fv;
           }, "Yield function")

      .def("df_ds",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_vec<double>(6);
            
            int ier = m.df_ds(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function gradient wrt. deviatoric stress")
      .def("df_dq",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_vec<double>(m.nhist());
            
            int ier = m.df_dq(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function gradient wrt. the history")

      .def("df_dsds",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_mat<double>(6,6);
            
            int ier = m.df_dsds(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function Hessian: stress-stress")

      .def("df_dsdq",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_mat<double>(6,m.nhist());
            
            int ier = m.df_dsdq(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function Hessian: stress-history")

      .def("df_dqds",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_mat<double>(m.nhist(),6);
            
            int ier = m.df_dqds(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function Hessian: history-stress")

      .def("df_dqdq",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_mat<double>(m.nhist(),m.nhist());
            
            int ier = m.df_dqdq(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function Hessian: history-history")

      .def("D",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto Dv = alloc_mat<double>(m.nhist(),m.nhist());
            
            int ier = m.D(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(Dv));

            return Dv;
           }, "Generalized plastic modulii")

      .def("D_inv",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto Dv = alloc_mat<double>(m.nhist(),m.nhist());
            
            int ier = m.D_inv(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(Dv));

            return Dv;
           }, "Inverse of generalized plastic modulii")
      ;
 
  py::class_<KinIsoJ2>(m, "KinIsoJ2", py::base<YieldSurface>())
      .def("K", &KinIsoJ2::K, "Isotropic hardening value")
      .def("Kp", &KinIsoJ2::Kp, "Isotropic hardening slope")
      .def("Hp", &KinIsoJ2::Hp, "Kinematic hardening slope")
      ;

  py::class_<LinearKinIsoJ2>(m, "LinearKinIsoJ2", py::base<KinIsoJ2>())
      .def(py::init<double, double, double>())
      .def_property_readonly("K0", &LinearKinIsoJ2::K0, "Initial yield stress")
      .def_property_readonly("Kb", &LinearKinIsoJ2::Kb, "Isotropic hardening slope")
      .def_property_readonly("Hb", &LinearKinIsoJ2::Hb, "Kinematic hardening slope")
      ;

  return m.ptr();
}


} // namespace neml
