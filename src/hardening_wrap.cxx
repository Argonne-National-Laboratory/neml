#include "hardening.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(hardening) {
  py::module m("hardening", "Various hardening rules.");

  py::class_<HardeningRule, std::shared_ptr<HardeningRule>>(m, "HardeningRule")
      .def_property_readonly("nhist", &HardeningRule::nhist, "Number of history variables.")

      .def("init_hist",
           [](const HardeningRule & m) -> py::array_t<double>
           {
            auto v = alloc_vec<double>(m.nhist());
            m.init_hist(arr2ptr<double>(v));
            return v;
           }, "Initialize history.")
        
      .def("q",
           [](const HardeningRule & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto q = alloc_vec<double>(m.nhist());
            m.q(arr2ptr<double>(alpha), T, arr2ptr<double>(q));
            return q;
           }, "Map alpha to q.")

      .def("dq_da",
           [](const HardeningRule & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto D = alloc_mat<double>(m.nhist(), m.nhist());
            m.dq_da(arr2ptr<double>(alpha), T, arr2ptr<double>(D));
            return D;
           }, "Gradient of map")
      ;
  
  py::class_<IsotropicHardeningRule, std::shared_ptr<IsotropicHardeningRule>>(m, "IsotropicHardeningRule", py::base<HardeningRule>())
      ;

  py::class_<LinearIsotropicHardeningRule, std::shared_ptr<LinearIsotropicHardeningRule>>(m, "LinearIsotropicHardeningRule", py::base<IsotropicHardeningRule>())
      .def(py::init<double, double>(), py::arg("s0"), py::arg("K"))
      .def(py::init<std::shared_ptr<Interpolate>, std::shared_ptr<Interpolate>>(), py::arg("s0"), py::arg("K"))
      
      .def("s0", &LinearIsotropicHardeningRule::s0)
      .def("K", &LinearIsotropicHardeningRule::K)
      ;

  py::class_<VoceIsotropicHardeningRule, std::shared_ptr<VoceIsotropicHardeningRule>>(m, "VoceIsotropicHardeningRule", py::base<IsotropicHardeningRule>())
      .def(py::init<double, double, double>(), py::arg("s0"), py::arg("R"), py::arg("d"))
      .def(py::init<std::shared_ptr<Interpolate>, std::shared_ptr<Interpolate>, std::shared_ptr<Interpolate>>(), py::arg("s0"), py::arg("R"), py::arg("d"))
      
      .def("s0", &VoceIsotropicHardeningRule::s0)
      .def("R", &VoceIsotropicHardeningRule::R)
      .def("d", &VoceIsotropicHardeningRule::d)
      ;

  py::class_<KinematicHardeningRule, std::shared_ptr<KinematicHardeningRule>>(m, "KinematicHardeningRule", py::base<HardeningRule>())
      ;

  py::class_<LinearKinematicHardeningRule, std::shared_ptr<LinearKinematicHardeningRule>>(m, "LinearKinematicHardeningRule", py::base<KinematicHardeningRule>())
      .def(py::init<double>(), py::arg("H"))
      .def(py::init<std::shared_ptr<Interpolate>>(), py::arg("H"))
      
      .def("H", &LinearKinematicHardeningRule::H)
      ;

  py::class_<CombinedHardeningRule, std::shared_ptr<CombinedHardeningRule>>(m, "CombinedHardeningRule", py::base<HardeningRule>())
      .def(py::init<std::shared_ptr<IsotropicHardeningRule>, std::shared_ptr<KinematicHardeningRule>>(),
           py::arg("iso"), py::arg("kin"))
      ;


  py::class_<NonAssociativeHardening, std::shared_ptr<NonAssociativeHardening>>(m, "NonAssociativeHardening")
      .def_property_readonly("ninter", &NonAssociativeHardening::ninter, "Number of q variables.")
      .def_property_readonly("nhist", &NonAssociativeHardening::nhist, "Number of a variables.")

      .def("init_hist",
           [](const NonAssociativeHardening & m) -> py::array_t<double>
           {
            auto v = alloc_vec<double>(m.nhist());
            int ier = m.init_hist(arr2ptr<double>(v));
            py_error(ier);
            return v;
           }, "Initialize history.")
        
      .def("q",
           [](const NonAssociativeHardening & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto q = alloc_vec<double>(m.ninter());
            int ier = m.q(arr2ptr<double>(alpha), T, arr2ptr<double>(q));
            py_error(ier);
            return q;
           }, "Map alpha to q.")

      .def("dq_da",
           [](const NonAssociativeHardening & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto D = alloc_mat<double>(m.ninter(), m.nhist());
            int ier = m.dq_da(arr2ptr<double>(alpha), T, arr2ptr<double>(D));
            py_error(ier);
            return D;
           }, "Gradient of map")

      .def("h",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            int ier = m.h(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule.")
      .def("dh_ds",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            int ier = m.dh_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule derivative with respect to stress.")
      .def("dh_da",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            int ier = m.dh_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule derivative with respect to history.")


      .def("h_time",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            int ier = m.h_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule, time part.")
      .def("dh_ds_time",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            int ier = m.dh_ds_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule (time) derivative with respect to stress.")
      .def("dh_da_time",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            int ier = m.dh_da_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule (time) derivative with respect to history.")

      .def("h_temp",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            int ier = m.h_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule, temperature part.")
      .def("dh_ds_temp",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            int ier = m.dh_ds_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule (temperature) derivative with respect to stress.")
      .def("dh_da_temp",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            int ier = m.dh_da_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Hardening rule (temperature) derivative with respect to history.")

      ;

      ;

  py::class_<GammaModel, std::shared_ptr<GammaModel>>(m, "GammaModel")
      .def("gamma", &GammaModel::gamma)
      .def("dgamma", &GammaModel::dgamma)
      ;

  py::class_<ConstantGamma, std::shared_ptr<ConstantGamma>>(m, "ConstantGamma", py::base<GammaModel>())
      .def(py::init<double>(), py::arg("g"))
      .def(py::init<std::shared_ptr<Interpolate>>(), py::arg("g"))

      .def("g", &ConstantGamma::g)
      ;

  py::class_<SatGamma, std::shared_ptr<SatGamma>>(m, "SatGamma", py::base<GammaModel>())
      .def(py::init<double,double,double>(), py::arg("gs"), py::arg("g0"),
           py::arg("beta"))
      .def(py::init<std::shared_ptr<Interpolate>,std::shared_ptr<Interpolate>,std::shared_ptr<Interpolate>>(), py::arg("gs"), py::arg("g0"),
           py::arg("beta"))
      .def("gs", &SatGamma::gs)
      .def("g0", &SatGamma::g0)
      .def("beta", &SatGamma::beta)
      ;

  py::class_<Chaboche, std::shared_ptr<Chaboche>>(m, "Chaboche", py::base<NonAssociativeHardening>())
      .def(py::init<std::shared_ptr<IsotropicHardeningRule>, std::vector<double>, std::vector<std::shared_ptr<GammaModel>>>(),
           py::arg("iso"), py::arg("c"), py::arg("gmodels"))
      .def(py::init<std::shared_ptr<IsotropicHardeningRule>, std::vector<std::shared_ptr<const Interpolate>>, std::vector<std::shared_ptr<GammaModel>>>(),
           py::arg("iso"), py::arg("c"), py::arg("gmodels"))

      .def(py::init<std::shared_ptr<IsotropicHardeningRule>, std::vector<double>, std::vector<std::shared_ptr<GammaModel>>, std::vector<double>, std::vector<double>>(),
           py::arg("iso"), py::arg("c"), py::arg("gmodels"), py::arg("A"), py::arg("a"))
      .def(py::init<std::shared_ptr<IsotropicHardeningRule>, std::vector<std::shared_ptr<const Interpolate>>, std::vector<std::shared_ptr<GammaModel>>, std::vector<std::shared_ptr<const Interpolate>>, std::vector<std::shared_ptr<const Interpolate>>>(),
           py::arg("iso"), py::arg("c"), py::arg("gmodels"), py::arg("A"), py::arg("a"))

      .def_property_readonly("n", &Chaboche::n, "Number of backstresses")
      .def("c",
         [](const Chaboche& m, double T) -> py::array_t<double>
         {
          auto cv = alloc_vec<double>(m.n());
          std::vector<double> vv = m.c(T);
          std::copy(vv.begin(), vv.end(), arr2ptr<double>(cv));
          return cv;
         }, "c material constant.")
      ;

  return m.ptr();

}


} // namespace neml
