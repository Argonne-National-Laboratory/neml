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
      
      .def_property_readonly("s0", &LinearIsotropicHardeningRule::s0)
      .def_property_readonly("K", &LinearIsotropicHardeningRule::K)
      ;

  py::class_<VoceIsotropicHardeningRule, std::shared_ptr<VoceIsotropicHardeningRule>>(m, "VoceIsotropicHardeningRule", py::base<IsotropicHardeningRule>())
      .def(py::init<double, double, double>(), py::arg("s0"), py::arg("R"), py::arg("d"))
      
      .def_property_readonly("s0", &VoceIsotropicHardeningRule::s0)
      .def_property_readonly("R", &VoceIsotropicHardeningRule::R)
      .def_property_readonly("d", &VoceIsotropicHardeningRule::d)
      ;

  py::class_<KinematicHardeningRule, std::shared_ptr<KinematicHardeningRule>>(m, "KinematicHardeningRule", py::base<HardeningRule>())
      ;

  py::class_<LinearKinematicHardeningRule, std::shared_ptr<LinearKinematicHardeningRule>>(m, "LinearKinematicHardeningRule", py::base<KinematicHardeningRule>())
      .def(py::init<double>(), py::arg("H"))
      
      .def_property_readonly("H", &LinearKinematicHardeningRule::H)
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
      ;

  py::class_<GammaModel, std::shared_ptr<GammaModel>>(m, "GammaModel")
      .def("gamma", &GammaModel::gamma)
      .def("dgamma", &GammaModel::dgamma)
      ;

  py::class_<ConstantGamma, std::shared_ptr<ConstantGamma>>(m, "ConstantGamma", py::base<GammaModel>())
      .def(py::init<double>(), py::arg("g"))
      .def_property_readonly("g", &ConstantGamma::g)
      ;

  py::class_<SatGamma, std::shared_ptr<SatGamma>>(m, "SatGamma", py::base<GammaModel>())
      .def(py::init<double,double,double>(), py::arg("gs"), py::arg("g0"),
           py::arg("beta"))
      .def_property_readonly("gs", &SatGamma::gs)
      .def_property_readonly("g0", &SatGamma::g0)
      .def_property_readonly("beta", &SatGamma::beta)
      ;

  py::class_<Chaboche, std::shared_ptr<Chaboche>>(m, "Chaboche", py::base<NonAssociativeHardening>())
      .def(py::init<std::shared_ptr<IsotropicHardeningRule>, std::vector<double>, std::vector<std::shared_ptr<GammaModel>>>(),
           py::arg("iso"), py::arg("c"), py::arg("gmodels"))
      .def("__init__",
           [](Chaboche & instance, std::shared_ptr<IsotropicHardeningRule> iso, py::array_t<double, py::array::c_style> c, py::array_t<double, py::array::c_style> r) 
           {
            if (c.request().ndim != 1) {
              throw std::runtime_error("c must be a vector!");
            }
            if (r.request().ndim != 1) {
              throw std::runtime_error("r must be a vector!");
            }
            if (c.request().shape[0] != r.request().shape[0]) {
              throw std::runtime_error("len(c) != len(r)!");
            }

            int n = c.request().shape[0];  

            new (&instance) Chaboche(iso, n, arr2ptr<double>(c), arr2ptr<double>(r));
           })
      .def_property_readonly("n", &Chaboche::n, "Number of backstresses")
      .def_property_readonly("c",
                             [](const Chaboche& m) -> py::array_t<double>
                             {
                              auto cv = alloc_vec<double>(m.n());
                              std::copy(m.c().begin(), m.c().end(), arr2ptr<double>(cv));
                              return cv;
                             }, "c material constant.")
      ;

  return m.ptr();

}


} // namespace neml
