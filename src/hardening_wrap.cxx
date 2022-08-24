#include "pyhelp.h" // include first to avoid annoying redef warning

#include "hardening.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(hardening, m) {
  py::module::import("neml.objects");

  m.doc() = "Various hardening rules.";

  py::class_<HardeningRule, HistoryNEMLObject, std::shared_ptr<HardeningRule>>(m, "HardeningRule")
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
  
  py::class_<IsotropicHardeningRule, HardeningRule, std::shared_ptr<IsotropicHardeningRule>>(m, "IsotropicHardeningRule")
      ;

  py::class_<LinearIsotropicHardeningRule, IsotropicHardeningRule, std::shared_ptr<LinearIsotropicHardeningRule>>(m, "LinearIsotropicHardeningRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<LinearIsotropicHardeningRule>(args,
                                                                         kwargs,
                                                                         {"s0",
                                                                         "K"});
                    }))      
      .def("s0", &LinearIsotropicHardeningRule::s0)
      .def("K", &LinearIsotropicHardeningRule::K)
      ;

  py::class_<InterpolatedIsotropicHardeningRule, IsotropicHardeningRule, std::shared_ptr<InterpolatedIsotropicHardeningRule>>(m, "InterpolatedIsotropicHardeningRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return
          create_object_python<InterpolatedIsotropicHardeningRule>(args,
                                                                   kwargs,
                                                                   {"flow"});
        }))
      ;

  py::class_<VoceIsotropicHardeningRule, IsotropicHardeningRule, std::shared_ptr<VoceIsotropicHardeningRule>>(m, "VoceIsotropicHardeningRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<VoceIsotropicHardeningRule>(args,
                                                                       kwargs,
                                                                       {"s0",
                                                                       "R", "d"});
                    }))

      .def("s0", &VoceIsotropicHardeningRule::s0)
      .def("R", &VoceIsotropicHardeningRule::R)
      .def("d", &VoceIsotropicHardeningRule::d)
      ;

  py::class_<PowerLawIsotropicHardeningRule, IsotropicHardeningRule, std::shared_ptr<PowerLawIsotropicHardeningRule>>(m, "PowerLawIsotropicHardeningRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<PowerLawIsotropicHardeningRule>(args,
                                                                       kwargs,
                                                                       {"s0",
                                                                       "A", "n"});
                    }))
      ;

  py::class_<CombinedIsotropicHardeningRule, IsotropicHardeningRule, std::shared_ptr<CombinedIsotropicHardeningRule>>(m, "CombinedIsotropicHardeningRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return
                      create_object_python<CombinedIsotropicHardeningRule>(args,
                                                                           kwargs,
                                                                           {"rules"});
                    }))
      .def_property_readonly("nrules", &CombinedIsotropicHardeningRule::nrules)
      ;

  py::class_<KinematicHardeningRule, HardeningRule, std::shared_ptr<KinematicHardeningRule>>(m, "KinematicHardeningRule")
      ;

  py::class_<LinearKinematicHardeningRule, KinematicHardeningRule, std::shared_ptr<LinearKinematicHardeningRule>>(m, "LinearKinematicHardeningRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<LinearKinematicHardeningRule>(args,
                                                                    kwargs, {"H"});
        }))
      .def("H", &LinearKinematicHardeningRule::H)
      ;

  py::class_<CombinedHardeningRule, HardeningRule, std::shared_ptr<CombinedHardeningRule>>(m, "CombinedHardeningRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<CombinedHardeningRule>(args, kwargs,
                                                             {"iso", "kin"});
        }))
      ;

  py::class_<NonAssociativeHardening, HistoryNEMLObject, std::shared_ptr<NonAssociativeHardening>>(m, "NonAssociativeHardening")
      .def_property_readonly("ninter", &NonAssociativeHardening::ninter, "Number of q variables.")
        
      .def("q",
           [](const NonAssociativeHardening & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto q = alloc_vec<double>(m.ninter());
            m.q(arr2ptr<double>(alpha), T, arr2ptr<double>(q));
            return q;
           }, "Map alpha to q.")

      .def("dq_da",
           [](const NonAssociativeHardening & m, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto D = alloc_mat<double>(m.ninter(), m.nhist());
            m.dq_da(arr2ptr<double>(alpha), T, arr2ptr<double>(D));
            return D;
           }, "Gradient of map")

      .def("h",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.h(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule.")
      .def("dh_ds",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.dh_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule derivative with respect to stress.")
      .def("dh_da",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            m.dh_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule derivative with respect to history.")


      .def("h_time",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.h_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule, time part.")
      .def("dh_ds_time",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.dh_ds_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule (time) derivative with respect to stress.")
      .def("dh_da_time",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            m.dh_da_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule (time) derivative with respect to history.")

      .def("h_temp",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.h_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule, temperature part.")
      .def("dh_ds_temp",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.dh_ds_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule (temperature) derivative with respect to stress.")
      .def("dh_da_temp",
           [](NonAssociativeHardening & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            m.dh_da_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule (temperature) derivative with respect to history.")
      ;

  py::class_<GammaModel, NEMLObject, std::shared_ptr<GammaModel>>(m, "GammaModel")
      .def("gamma", &GammaModel::gamma)
      .def("dgamma", &GammaModel::dgamma)
      ;

  py::class_<ConstantGamma, GammaModel, std::shared_ptr<ConstantGamma>>(m, "ConstantGamma")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ConstantGamma>(args, kwargs, {"g"});
        }))
      .def("g", &ConstantGamma::g)
      ;

  py::class_<SatGamma, GammaModel, std::shared_ptr<SatGamma>>(m, "SatGamma")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<SatGamma>(args, kwargs, {"gs", "g0",
                                                "beta"});
        }))
      .def("gs", &SatGamma::gs)
      .def("g0", &SatGamma::g0)
      .def("beta", &SatGamma::beta)
      ;

  py::class_<Chaboche, NonAssociativeHardening, std::shared_ptr<Chaboche>>(m, "Chaboche")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<Chaboche>(args, kwargs, {"iso", "C",
                                                "gmodels", "A", "a"});
        }))
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

  py::class_<ChabocheVoceRecovery, NonAssociativeHardening,
      std::shared_ptr<ChabocheVoceRecovery>>(m, "ChabocheVoceRecovery")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ChabocheVoceRecovery>(args, kwargs, {"s0",
                                                            "theta0", "Rmax",
                                                            "Rmin", "r1", "r2", "C",
                                                "gmodels", "A", "a"});
        }))
      .def_property_readonly("n", &ChabocheVoceRecovery::n, "Number of backstresses")
      ;
}

} // namespace neml
