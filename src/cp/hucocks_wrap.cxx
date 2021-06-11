#include "../pyhelp.h"

#include "hucocks.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(hucocks, m) {
  py::module::import("neml.objects");

  m.doc() = "Objects for the Hu & Cocks 316H model";

  py::class_<HuCocksPrecipitationModel, NEMLObject,
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
    .def("populate_history", &HuCocksPrecipitationModel::populate_history)
    .def("init_history", &HuCocksPrecipitationModel::init_history)
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
    ;

} // PYBIND!!_MODULE


}
