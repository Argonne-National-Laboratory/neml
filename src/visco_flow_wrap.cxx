#include "pyhelp.h" // include first to avoid annoying redef warning

#include "visco_flow.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {
PYBIND11_MODULE(visco_flow, m) {
  py::module::import("neml.objects");

  m.doc() = "Viscoplastic flow models.";

  py::class_<ViscoPlasticFlowRule, HistoryNEMLObject, std::shared_ptr<ViscoPlasticFlowRule>>(m, "ViscoPlasticFlowRule")
      .def("y",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> double
           {
            double fv;
            m.y(arr2ptr<double>(s), arr2ptr<double>(alpha), T, fv);
            return fv;
           }, "Plastic multiplier value.")
      .def("dy_ds",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            m.dy_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Plastic multiplier derivative with respect to stress.")
      .def("dy_da",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.dy_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Plastic multiplier derivative with respect to history.")

      .def("g",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            m.g(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule, rate part.")
      .def("dg_ds",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            m.dg_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule (rate) derivative with respect to stress.")
      .def("dg_da",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,m.nhist());
            m.dg_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule (rate) derivative with respect to history.")

      .def("g_time",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            m.g_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule, time part.")
      .def("dg_ds_time",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            m.dg_ds_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule (time) derivative with respect to stress.")
      .def("dg_da_time",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,m.nhist());
            m.dg_da_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule (time) derivative with respect to history.")

      .def("g_temp",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            m.g_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule, temperature part.")
      .def("dg_ds_temp",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,6);
            m.dg_ds_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule (temperature) derivative with respect to stress.")
      .def("dg_da_temp",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(6,m.nhist());
            m.dg_da_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Flow rule (temperature) derivative with respect to history.")

      .def("h",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.h(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule.")
      .def("dh_ds",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.dh_ds(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule derivative with respect to stress.")
      .def("dh_da",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            m.dh_da(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule derivative with respect to history.")

      .def("h_time",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.h_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule, time part.")
      .def("dh_ds_time",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.dh_ds_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule (time) derivative with respect to stress.")
      .def("dh_da_time",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            m.dh_da_time(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule (time) derivative with respect to history.")

      .def("h_temp",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(m.nhist());
            m.h_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule, temperature part.")
      .def("dh_ds_temp",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),6);
            m.dh_ds_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule (temperature) derivative with respect to stress.")
      .def("dh_da_temp",
           [](ViscoPlasticFlowRule & m, py::array_t<double, py::array::c_style> s, py::array_t<double, py::array::c_style> alpha, double T) -> py::array_t<double>
           {
            auto f = alloc_mat<double>(m.nhist(),m.nhist());
            m.dh_da_temp(arr2ptr<double>(s), arr2ptr<double>(alpha), T, arr2ptr<double>(f));
            return f;
           }, "Hardening rule (temperature) derivative with respect to history.")

      ;

  py::class_<SuperimposedViscoPlasticFlowRule, ViscoPlasticFlowRule, std::shared_ptr<SuperimposedViscoPlasticFlowRule>>(m, "SuperimposedViscoPlasticFlowRule")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<SuperimposedViscoPlasticFlowRule>(args, kwargs, {"flow_rules"});
      }))
      ;

  py::class_<LinearViscousFlow, ViscoPlasticFlowRule, std::shared_ptr<LinearViscousFlow>>(m, "LinearViscousFlow")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<LinearViscousFlow>(args, kwargs, {"surface",
                                                     "eta"});
      }))
      ;

  py::class_<GFlow, NEMLObject, std::shared_ptr<GFlow>>(m, "GFlow")
      .def("g", &GFlow::g, "g function in Perzyna model")
      .def("dg", &GFlow::dg, "Derivative of g wrt f")
      ;

  py::class_<GPowerLaw, GFlow, std::shared_ptr<GPowerLaw>>(m, "GPowerLaw")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<GPowerLaw>(args, kwargs, {"n", "eta"});
        }))

      .def("n", &GPowerLaw::n)
      .def("eta", &GPowerLaw::eta)
      ;

  py::class_<PerzynaFlowRule, ViscoPlasticFlowRule, std::shared_ptr<PerzynaFlowRule>>(m, "PerzynaFlowRule")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<PerzynaFlowRule>(args, kwargs, {"surface",
                                                     "hardening",
                                                     "g"});
      }))
      ;

  py::class_<FluidityModel, NEMLObject, std::shared_ptr<FluidityModel>>(m, "FluidityModel")
      .def("eta", &FluidityModel::eta)
      .def("deta", &FluidityModel::deta)
      ;

  py::class_<ConstantFluidity, FluidityModel, std::shared_ptr<ConstantFluidity>>(m, "ConstantFluidity")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ConstantFluidity>(args, kwargs, {"eta"});
        }))
      ;

  py::class_<SaturatingFluidity, FluidityModel, std::shared_ptr<SaturatingFluidity>>(m, "SaturatingFluidity")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<SaturatingFluidity>(args, kwargs, {"K0", "A", "b"});
        }))
      ;

  py::class_<ChabocheFlowRule, ViscoPlasticFlowRule, std::shared_ptr<ChabocheFlowRule>>(m, "ChabocheFlowRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<ChabocheFlowRule>(args, kwargs, {"surface", "hardening",
                                                        "fluidity", "n"});
        }))
      ;

  py::class_<YaguchiGr91FlowRule, ViscoPlasticFlowRule, std::shared_ptr<YaguchiGr91FlowRule>>(m, "YaguchiGr91FlowRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<YaguchiGr91FlowRule>(args, kwargs, {});
        }))
      .def("D", &YaguchiGr91FlowRule::D)
      .def("n", &YaguchiGr91FlowRule::n)
      .def("a10", &YaguchiGr91FlowRule::a10)
      .def("C2", &YaguchiGr91FlowRule::C2)
      .def("a2", &YaguchiGr91FlowRule::a2)
      .def("g1", &YaguchiGr91FlowRule::g1)
      .def("g2", &YaguchiGr91FlowRule::g2)
      .def("m", &YaguchiGr91FlowRule::m)
      .def("br", &YaguchiGr91FlowRule::br)
      .def("bh", &YaguchiGr91FlowRule::bh)
      .def("A", &YaguchiGr91FlowRule::A)
      .def("B", &YaguchiGr91FlowRule::B)
      .def("d", &YaguchiGr91FlowRule::d)
      .def("q", &YaguchiGr91FlowRule::q)
      .def("C1", &YaguchiGr91FlowRule::C1)
      ;
}

} // namespace neml
