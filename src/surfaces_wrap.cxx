#include "surfaces.h"

#include "pyhelp.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(surfaces) {
  py::module m("surfaces", "Various yield surfaces.");

  py::class_<YieldSurface, std::shared_ptr<YieldSurface>>(m, "YieldSurface")
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
      ;
 
  py::class_<IsoJ2, std::shared_ptr<IsoJ2>>(m, "IsoJ2", py::base<YieldSurface>())
      .def(py::init<>())
      ;

  return m.ptr();
}


} // namespace neml
