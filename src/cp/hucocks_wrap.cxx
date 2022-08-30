#include "pyhelp.h"

#include "cp/hucocks.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(hucocks, m) {
  py::module::import("neml.objects");
  py::module::import("neml.cp.slipharden");
  py::module::import("neml.cp.sliprules");

  m.doc() = "Objects for the Hu & Cocks 316H model";

  py::class_<HuCocksPrecipitationModel, HistoryNEMLObject,
      std::shared_ptr<HuCocksPrecipitationModel>>(m,
                                                  "HuCocksPrecipitationModel")
    .def(py::init([](py::args args, py::kwargs kwargs)
                  {
                    return create_object_python<HuCocksPrecipitationModel>(
                        args, kwargs, {"c0", "cp", "ceq", "am", "N0", "Vm",
                        "chi", "D0", "Q0", "Cf"});
                  }))
    .def("varnames", &HuCocksPrecipitationModel::varnames)
    .def("set_varnames", &HuCocksPrecipitationModel::set_varnames)
    .def("f", &HuCocksPrecipitationModel::f)
    .def("r", &HuCocksPrecipitationModel::r)
    .def("N", &HuCocksPrecipitationModel::N)
    .def("rate", &HuCocksPrecipitationModel::rate)
    .def("jac", &HuCocksPrecipitationModel::jac)
    .def("f_rate", &HuCocksPrecipitationModel::f_rate)
    .def("df_df", &HuCocksPrecipitationModel::df_df)
    .def("df_dr", &HuCocksPrecipitationModel::df_dr)
    .def("df_dN", &HuCocksPrecipitationModel::df_dN)
    .def("r_rate", &HuCocksPrecipitationModel::r_rate)
    .def("dr_df", &HuCocksPrecipitationModel::dr_df)
    .def("dr_dr", &HuCocksPrecipitationModel::dr_dr)
    .def("dr_dN", &HuCocksPrecipitationModel::dr_dN)
    .def("N_rate", &HuCocksPrecipitationModel::N_rate)
    .def("dN_df", &HuCocksPrecipitationModel::dN_df)
    .def("dN_dr", &HuCocksPrecipitationModel::dN_dr)
    .def("dN_dN", &HuCocksPrecipitationModel::dN_dN)
    .def_property_readonly("nspecies", &HuCocksPrecipitationModel::nspecies)
    .def("c", &HuCocksPrecipitationModel::c)
    .def("dc_df", &HuCocksPrecipitationModel::dc_df)
    .def("Gv", &HuCocksPrecipitationModel::Gv)
    .def("dG_df", &HuCocksPrecipitationModel::dG_df)
    .def_property_readonly("vm", &HuCocksPrecipitationModel::vm)
    .def_property_readonly("fs", &HuCocksPrecipitationModel::fs)
    .def_property_readonly("rs", &HuCocksPrecipitationModel::rs)
    .def_property_readonly("Ns", &HuCocksPrecipitationModel::Ns)
    ;

  py::class_<DislocationSpacingHardening, SlipHardening, std::shared_ptr<DislocationSpacingHardening>>(m, "DislocationSpacingHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<DislocationSpacingHardening>(
                          args, kwargs, {"J1", "J2", "K", "L0", "a", "b", "G",
                          "L"});
                    }))
      ;

  py::class_<HuCocksHardening, SlipHardening, std::shared_ptr<HuCocksHardening>>(m, "HuCocksHardening")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<HuCocksHardening>(
                          args, kwargs, {"dmodel", "pmodels", "ap", "ac", "b",
                          "G"});
                    }))
      ;
  py::class_<ArrheniusSlipRule, SlipStrengthSlipRule,
        std::shared_ptr<ArrheniusSlipRule>>(m, "ArrheniusSlipRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<ArrheniusSlipRule>(
                          args, kwargs, {"resistance", "g0", "A",
                          "B", "b", "a0", "G0"});
                    }))
      ;

} // PYBIND!!_MODULE


}
