#include "../pyhelp.h"

#include "sliprules.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(sliprules, m) {
  m.doc() = "Crystal plasticity slip rate relations";

  py::class_<SlipRule, NEMLObject, std::shared_ptr<SlipRule>>(m, "SlipRule")
      .def("populate_history", &SlipRule::populate_history)
      .def("init_history", &SlipRule::init_history)
      .def("slip", &SlipRule::slip)
      .def("d_slip_d_s", &SlipRule::d_slip_d_s)
      .def("d_slip_d_h", &SlipRule::d_slip_d_h)
      .def("hist_rate", &SlipRule::hist_rate)
      .def("d_hist_rate_d_stress", &SlipRule::d_hist_rate_d_stress)
      .def("d_hist_rate_d_hist", &SlipRule::d_hist_rate_d_hist)
      ;

  py::class_<SlipStrengthSlipRule, SlipRule,
        std::shared_ptr<SlipStrengthSlipRule>>(m, "SlipStrengthSlipRule")
      .def("sslip", &SlipStrengthSlipRule::sslip)
      .def("d_sslip_dtau", &SlipStrengthSlipRule::d_sslip_dtau)
      .def("d_sslip_dstrength", &SlipStrengthSlipRule::d_sslip_dstrength)
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
