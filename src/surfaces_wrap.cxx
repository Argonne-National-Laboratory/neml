#include "pyhelp.h" // include first to avoid annoying redef warning

#include "surfaces.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(surfaces, m) {
  py::module::import("neml.objects");

  m.doc() = "Yield (and flow) surface definitions.";

  py::class_<YieldSurface, NEMLObject, std::shared_ptr<YieldSurface>>(m, "YieldSurface")
      .def_property_readonly("nhist", &YieldSurface::nhist, "Number of history variables.")

      .def("f", 
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> double
           {
            double fv;

            m.f(arr2ptr<double>(s), arr2ptr<double>(h), T, fv);

            return fv;
           }, "Yield function")

      .def("df_ds",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_vec<double>(6);
            
            m.df_ds(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function gradient wrt. deviatoric stress")
      .def("df_dq",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_vec<double>(m.nhist());
            
            m.df_dq(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function gradient wrt. the history")

      .def("df_dsds",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_mat<double>(6,6);
            
            m.df_dsds(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function Hessian: stress-stress")

      .def("df_dsdq",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_mat<double>(6,m.nhist());
            
            m.df_dsdq(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function Hessian: stress-history")

      .def("df_dqds",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_mat<double>(m.nhist(),6);
            
            m.df_dqds(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function Hessian: history-stress")

      .def("df_dqdq",
           [](const YieldSurface & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> h, double T) -> py::array_t<double>
           {
            auto deriv = alloc_mat<double>(m.nhist(),m.nhist());
            
            m.df_dqdq(arr2ptr<double>(s), arr2ptr<double>(h), T, arr2ptr<double>(deriv));

            return deriv;
           }, "Yield function Hessian: history-history")
      ;
 
  py::class_<IsoJ2, YieldSurface, std::shared_ptr<IsoJ2>>(m, "IsoJ2")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<IsoJ2>(args, kwargs, {});
                    }))
      ;

  py::class_<IsoKinJ2, YieldSurface, std::shared_ptr<IsoKinJ2>>(m, "IsoKinJ2")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<IsoKinJ2>(args, kwargs, {});
                    }))
      ;

  py::class_<IsoKinJ2I1, YieldSurface, std::shared_ptr<IsoKinJ2I1>>(m, "IsoKinJ2I1")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<IsoKinJ2I1>(args, kwargs, 
                                                              {"h", "l"});
                    }))
      ;

  py::class_<IsoJ2I1, YieldSurface, std::shared_ptr<IsoJ2I1>>(m, "IsoJ2I1")
       .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<IsoJ2I1>(args, kwargs, 
                                                              {"h", "l"});
                    })) 
      ;
}

} // namespace neml
