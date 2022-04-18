#include "pyhelp.h" // include first to avoid annoying redef warning

#include "elasticity.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(elasticity, m) {
  py::module::import("neml.objects");

  m.doc() = "Elastic models.";

  py::class_<LinearElasticModel, NEMLObject, std::shared_ptr<LinearElasticModel>>(m, "LinearElasticModel")
      .def("C",
           [](const LinearElasticModel & m, double T) -> py::array_t<double>
           {
            auto C = alloc_mat<double>(6,6);
            m.C(T, arr2ptr<double>(C));
            return C;
           }, "Return stiffness elasticity matrix.")

      .def("S",
           [](const LinearElasticModel & m, double T) -> py::array_t<double>
           {
            auto S = alloc_mat<double>(6,6);
            m.S(T, arr2ptr<double>(S));
            return S;
           }, "Return compliance elasticity matrix.")
      .def("C_tensor",
           [](const LinearElasticModel & m, double T) -> SymSymR4
           {
            return m.C(T);
           }, "Return stiffness elasticity tensor.")

      .def("S_tensor",
           [](const LinearElasticModel & m, double T) -> SymSymR4
           {
            return m.S(T);
           }, "Return compliance elasticity tensor.")
      .def("C_tensor",
           [](const LinearElasticModel & m, double T, const Orientation & Q) -> SymSymR4
           {
            return m.C(T, Q);
           }, "Return rotated stiffness elasticity tensor.")

      .def("S_tensor",
           [](const LinearElasticModel & m, double T, const Orientation & Q) -> SymSymR4
           {
            return m.S(T, Q);
           }, "Return rotated compliance elasticity tensor.")
      .def("G", (double (LinearElasticModel::*)(double) const) &LinearElasticModel::G)
      .def("G", (double (LinearElasticModel::*)(double, const Orientation &,
                                                const Vector &, const Vector &)
                 const)
           &LinearElasticModel::G)
      ;

  py::class_<IsotropicLinearElasticModel, LinearElasticModel, std::shared_ptr<IsotropicLinearElasticModel>>(m, "IsotropicLinearElasticModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<IsotropicLinearElasticModel>(args, kwargs, {"m1", "m1_type", "m2", "m2_type"});
        }))
      .def("E", &IsotropicLinearElasticModel::E, "Young's modulus as a function of temperature.")
      .def("nu", &IsotropicLinearElasticModel::nu, "Poisson's ratio as a function of temperature.")
      .def("K", &IsotropicLinearElasticModel::K, "Bulk modulus as a function of temperature.")
      ;

  py::class_<CubicLinearElasticModel, LinearElasticModel, std::shared_ptr<CubicLinearElasticModel>>(m, "CubicLinearElasticModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<CubicLinearElasticModel>(args, kwargs, {"m1", "m2", "m3", "method"});
        }))
  ;

  py::class_<TransverseIsotropicLinearElasticModel, LinearElasticModel,
      std::shared_ptr<TransverseIsotropicLinearElasticModel>>(m,
                                                              "TransverseIsotropicLinearElasticModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return
          create_object_python<TransverseIsotropicLinearElasticModel>(args,
                                                                      kwargs,
                                                                      {"m1",
                                                                      "m2",
                                                                      "m3",
                                                                      "m4", "m5",  "method"});
        }))
  ;
}

}
