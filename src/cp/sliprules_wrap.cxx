#include "pyhelp.h"

#include "cp/sliprules.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(sliprules, m) {
  m.doc() = "Crystal plasticity slip rate relations";

  py::class_<SlipRule, HistoryNEMLObject, std::shared_ptr<SlipRule>>(m, "SlipRule")
      .def("strength", &SlipRule::strength)
      .def("slip", &SlipRule::slip)
      .def("d_slip_d_s", &SlipRule::d_slip_d_s)
      .def("d_slip_d_h", &SlipRule::d_slip_d_h)
      .def("hist_rate", &SlipRule::hist_rate)
      .def("d_hist_rate_d_stress", &SlipRule::d_hist_rate_d_stress)
      .def("d_hist_rate_d_hist", &SlipRule::d_hist_rate_d_hist)
      .def("sum_slip", &SlipRule::sum_slip)
      .def("d_sum_slip_d_stress", &SlipRule::d_sum_slip_d_stress)
      .def("d_sum_slip_d_hist", &SlipRule::d_sum_slip_d_hist)
      .def_property_readonly("use_nye", &SlipRule::use_nye)
      ;

  py::class_<SlipMultiStrengthSlipRule, SlipRule,
        std::shared_ptr<SlipMultiStrengthSlipRule>>(m, "SlipMultiStrengthSlipRule")
      .def("sslip", &SlipMultiStrengthSlipRule::sslip)
      .def("d_sslip_dtau", &SlipMultiStrengthSlipRule::d_sslip_dtau)
      .def("d_sslip_dstrength", &SlipMultiStrengthSlipRule::d_sslip_dstrength)
      ;

  py::class_<KinematicPowerLawSlipRule, SlipMultiStrengthSlipRule,
        std::shared_ptr<KinematicPowerLawSlipRule>>(m, "KinematicPowerLawSlipRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<KinematicPowerLawSlipRule>(
                          args, kwargs, {"backstrength", "isostrength",
                          "flowresistance", "gamma0", "n"});
                    }))
      ;

  py::class_<SlipStrengthSlipRule, SlipMultiStrengthSlipRule,
        std::shared_ptr<SlipStrengthSlipRule>>(m, "SlipStrengthSlipRule")
      .def("scalar_sslip", &SlipStrengthSlipRule::scalar_sslip)
      .def("scalar_d_sslip_dtau", &SlipStrengthSlipRule::scalar_d_sslip_dtau)
      .def("scalar_d_sslip_dstrength",
           &SlipStrengthSlipRule::scalar_d_sslip_dstrength)
      ;

  py::class_<PowerLawSlipRule, SlipStrengthSlipRule,
        std::shared_ptr<PowerLawSlipRule>>(m, "PowerLawSlipRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<PowerLawSlipRule>(
                          args, kwargs, {"resistance", "gamma0", "n"});
                    }))
      ;

} // PYBIND11_MODULE

} // namespace neml
