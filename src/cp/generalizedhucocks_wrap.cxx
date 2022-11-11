#include "pyhelp.h"

#include "cp/generalizedhucocks.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml
{

PYBIND11_MODULE(generalizedhucocks, m)
{
  py::module::import("neml.objects");
  py::module::import("neml.cp.slipharden");

  m.doc() = "Objects for the Generalized Hu & Cocks 316H model";

  py::class_<GeneralizedHuCocksSpecies, NEMLObject, std::shared_ptr<GeneralizedHuCocksSpecies>>(
      m, "GeneralizedHuCocksSpecies")
      .def(py::init(
          [](py::args args, py::kwargs kwargs)
          {
            return create_object_python<GeneralizedHuCocksSpecies>(
                args, kwargs, {"composition", "c0", "ceq"});
          }))
      .def_readonly("composition", &GeneralizedHuCocksSpecies::composition)
      .def_readonly("c0", &GeneralizedHuCocksSpecies::c0)
      .def_readonly("ceq", &GeneralizedHuCocksSpecies::ceq);

  py::class_<GeneralizedHuCocksPrecipitate,
             HistoryNEMLObject,
             std::shared_ptr<GeneralizedHuCocksPrecipitate>>(m, "GeneralizedHuCocksPrecipitate")
      .def(py::init(
          [](py::args args, py::kwargs kwargs)
          {
            return create_object_python<GeneralizedHuCocksPrecipitate>(args,
                                                                       kwargs,
                                                                       {"composition",
                                                                        "species",
                                                                        "rate",
                                                                        "cp",
                                                                        "am",
                                                                        "Vm",
                                                                        "D0",
                                                                        "Q0",
                                                                        "N0",
                                                                        "chi",
                                                                        "Cf"});
          }))
      .def("f", &GeneralizedHuCocksPrecipitate::f)
      .def("d_f_d_r", &GeneralizedHuCocksPrecipitate::d_f_d_r)
      .def("d_f_d_N", &GeneralizedHuCocksPrecipitate::d_f_d_N)
      .def("r", &GeneralizedHuCocksPrecipitate::r)
      .def("N", &GeneralizedHuCocksPrecipitate::N)
      .def("rs", &GeneralizedHuCocksPrecipitate::rs)
      .def("Ns", &GeneralizedHuCocksPrecipitate::Ns)
      .def_readonly("composition", &GeneralizedHuCocksPrecipitate::composition)
      .def_readonly("species", &GeneralizedHuCocksPrecipitate::species)
      .def_readonly("species_names", &GeneralizedHuCocksPrecipitate::species_names)
      .def_readonly("cp", &GeneralizedHuCocksPrecipitate::cp)
      .def_readonly("am", &GeneralizedHuCocksPrecipitate::am)
      .def_readonly("Vm", &GeneralizedHuCocksPrecipitate::Vm)
      .def_readonly("D0", &GeneralizedHuCocksPrecipitate::D0)
      .def_readonly("Q0", &GeneralizedHuCocksPrecipitate::Q0)
      .def_readonly("N0", &GeneralizedHuCocksPrecipitate::N0)
      .def_readonly("chi", &GeneralizedHuCocksPrecipitate::chi)
      .def_readonly("Cf", &GeneralizedHuCocksPrecipitate::Cf);

  py::class_<GeneralizedHuCocksPrecipitationModel,
             HistoryNEMLObject,
             std::shared_ptr<GeneralizedHuCocksPrecipitationModel>>(
      m, "GeneralizedHuCocksPrecipitationModel")
      .def(py::init(
          [](py::args args, py::kwargs kwargs)
          {
            return create_object_python<GeneralizedHuCocksPrecipitationModel>(
                args, kwargs, {"species", "precipitates"});
          }))
      .def("diffusivity", &GeneralizedHuCocksPrecipitationModel::diffusivity)
      .def("concentration", &GeneralizedHuCocksPrecipitationModel::concentration)
      .def("d_concentration_d_f", &GeneralizedHuCocksPrecipitationModel::d_concentration_d_f)
      .def("gibbs_free_energy", &GeneralizedHuCocksPrecipitationModel::gibbs_free_energy)
      .def("d_gibbs_free_energy_d_f",
           &GeneralizedHuCocksPrecipitationModel::d_gibbs_free_energy_d_f)
      .def("growth_rate", &GeneralizedHuCocksPrecipitationModel::growth_rate)
      .def("d_growth_rate", &GeneralizedHuCocksPrecipitationModel::d_growth_rate)
      .def("ripening_rate", &GeneralizedHuCocksPrecipitationModel::ripening_rate)
      .def("d_ripening_rate", &GeneralizedHuCocksPrecipitationModel::d_ripening_rate)
      .def("switching_function", &GeneralizedHuCocksPrecipitationModel::switching_function)
      .def("mixed_rate", &GeneralizedHuCocksPrecipitationModel::mixed_rate)
      .def("d_mixed_rate", &GeneralizedHuCocksPrecipitationModel::d_mixed_rate)
      .def("rate", &GeneralizedHuCocksPrecipitationModel::rate)
      .def("d_rate", &GeneralizedHuCocksPrecipitationModel::d_rate)
      .def_readonly("species", &GeneralizedHuCocksPrecipitationModel::species)
      .def_readonly("precipitates", &GeneralizedHuCocksPrecipitationModel::precipitates)
      .def_readonly("kboltz", &GeneralizedHuCocksPrecipitationModel::kboltz)
      .def_readonly("Na", &GeneralizedHuCocksPrecipitationModel::Na)
      .def_readonly("R", &GeneralizedHuCocksPrecipitationModel::R);

  py::class_<GeneralizedHuCocksHardening,
             SlipHardening,
             std::shared_ptr<GeneralizedHuCocksHardening>>(m, "GeneralizedHuCocksHardening")
      .def(py::init(
          [](py::args args, py::kwargs kwargs)
          {
            return create_object_python<GeneralizedHuCocksHardening>(
                args, kwargs, {"dmodel", "pmodel", "ap", "ac", "b", "G"});
          }))
      .def("dmodel", &GeneralizedHuCocksHardening::dmodel)
      .def("pmodel", &GeneralizedHuCocksHardening::pmodel)
      .def("hist_to_tau", &GeneralizedHuCocksHardening::hist_to_tau)
      .def("d_hist_to_tau", &GeneralizedHuCocksHardening::d_hist_to_tau);
} // PYBIND!!_MODULE

}
