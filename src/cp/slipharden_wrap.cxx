#include "pyhelp.h"

#include "cp/slipharden.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(slipharden, m) {
  m.doc() = "Crystal plasticity slip rate relations";

  py::class_<SlipHardening, HistoryNEMLObject, std::shared_ptr<SlipHardening>>(m, "SlipHardening")
      .def_property_readonly("varnames", &SlipHardening::varnames)
      .def("set_varnames", &SlipHardening::set_varnames)
      .def("hist_to_tau", &SlipHardening::hist_to_tau)
      .def("d_hist_to_tau", &SlipHardening::d_hist_to_tau)
      .def("hist", &SlipHardening::hist)
      .def("d_hist_d_s", &SlipHardening::d_hist_d_s)
      .def("d_hist_d_h", &SlipHardening::d_hist_d_h)
      .def("d_hist_d_h_ext", &SlipHardening::d_hist_d_h_ext)
      .def_property_readonly("use_nye", &SlipHardening::use_nye)
      ;

  py::class_<FixedStrengthHardening, SlipHardening, std::shared_ptr<FixedStrengthHardening>>(m, "FixedStrengthHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<FixedStrengthHardening>(
                          args, kwargs, {"strengths"});
                    }))
    ;

  py::class_<VocePerSystemHardening, SlipHardening,
      std::shared_ptr<VocePerSystemHardening>>(m, "VocePerSystemHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<VocePerSystemHardening>(
                          args, kwargs, {"initial", "k", "saturation",
                          "m"});
                    }))
    ;

  py::class_<FASlipHardening, SlipHardening,
      std::shared_ptr<FASlipHardening>>(m, "FASlipHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<FASlipHardening>(
                          args, kwargs, {"k", "saturation"});
                    }))
    ;

  py::class_<GeneralLinearHardening, SlipHardening, std::shared_ptr<GeneralLinearHardening>>(m, "GeneralLinearHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<GeneralLinearHardening>(
                          args, kwargs, {"M", "tau_0"});
                    }))
      ;

  py::class_<SimpleLinearHardening, SlipHardening,
      std::shared_ptr<SimpleLinearHardening>>(m, "SimpleLinearHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<SimpleLinearHardening>(
                          args, kwargs, {"G", "tau_0"});
                    }))
      ;

  py::class_<LANLTiModel, SlipHardening,
      std::shared_ptr<LANLTiModel>>(m, "LANLTiModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<LANLTiModel>(
                          args, kwargs, {"tau_0", "C_st", "mu", 
						  "k1", "k2"});
                    }))
      ;


  py::class_<SlipSingleHardening, SlipHardening,
        std::shared_ptr<SlipSingleHardening>>(m, "SlipSingleHardening")
      .def("hist_map", &SlipSingleHardening::hist_map)
      .def("d_hist_map", &SlipSingleHardening::d_hist_map)
      ;

  py::class_<SumSlipSingleStrengthHardening, SlipSingleHardening,
      std::shared_ptr<SumSlipSingleStrengthHardening>>(m,
                                                       "SumSlipSingleStrengthHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<SumSlipSingleStrengthHardening>(
                          args, kwargs, {"models"});
                    }))
      ;

  py::class_<SlipSingleStrengthHardening, SlipSingleHardening,
      std::shared_ptr<SlipSingleStrengthHardening>>(m,
                                                    "SlipSingleStrengthHardening")
      .def("init_strength", &SlipSingleStrengthHardening::init_strength)
      .def("hist_rate", &SlipSingleStrengthHardening::hist_rate)
      .def("d_hist_rate_d_stress",
           &SlipSingleStrengthHardening::d_hist_rate_d_stress)
      .def("d_hist_rate_d_hist",
           &SlipSingleStrengthHardening::d_hist_rate_d_hist)
      .def("d_hist_rate_d_hist_ext",
           &SlipSingleStrengthHardening::d_hist_rate_d_hist_ext)
      .def("static_strength", &SlipSingleStrengthHardening::static_strength)
      .def("nye_contribution", &SlipSingleStrengthHardening::nye_contribution)
      .def("nye_part", &SlipSingleStrengthHardening::nye_part)
      ;

  py::class_<PlasticSlipHardening, SlipSingleStrengthHardening,
        std::shared_ptr<PlasticSlipHardening>>(m, "PlasticSlipHardening")
      .def("hist_factor", &PlasticSlipHardening::hist_factor)
      .def("d_hist_factor", &PlasticSlipHardening::d_hist_factor)
      ;

  py::class_<VoceSlipHardening, PlasticSlipHardening,
        std::shared_ptr<VoceSlipHardening>>(m, "VoceSlipHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<VoceSlipHardening>(
                          args, kwargs, {"tau_sat", "b", "tau_0"});
                    }))
      ;

  py::class_<LinearSlipHardening, PlasticSlipHardening,
        std::shared_ptr<LinearSlipHardening>>(m, "LinearSlipHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<LinearSlipHardening>(
                          args, kwargs, {"tau0", "k1", "k2"});
                    }))
      ;

} // PYBIND11_MODULE

} // namespace neml
