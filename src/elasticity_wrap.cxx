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
            int ier = m.C(T, arr2ptr<double>(C));
            py_error(ier);
            return C;
           }, "Return stiffness elasticity matrix.")

      .def("S",
           [](const LinearElasticModel & m, double T) -> py::array_t<double>
           {
            auto S = alloc_mat<double>(6,6);
            int ier = m.S(T, arr2ptr<double>(S));
            py_error(ier);
            return S;
           }, "Return compliance elasticity matrix.")
      .def("E", &LinearElasticModel::E, "Young's modulus as a function of temperature.")
      .def("nu", &LinearElasticModel::nu, "Poisson's ratio as a function of temperature.")
      .def("G", &LinearElasticModel::G, "Shear modulus as a function of temperature.")
      .def("K", &LinearElasticModel::K, "Bulk modulus as a function of temperature.")
      .def_property_readonly("valid", &LinearElasticModel::valid, "Good or dummy model.")
      ;

  py::class_<IsotropicLinearElasticModel, LinearElasticModel, std::shared_ptr<IsotropicLinearElasticModel>>(m, "IsotropicLinearElasticModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<IsotropicLinearElasticModel>(args, kwargs, {"m1", "m1_type", "m2", "m2_type"});
        }))
      ;
  
  py::class_<BlankElasticModel, LinearElasticModel, std::shared_ptr<BlankElasticModel>>(m, "BlankElasticModel")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<BlankElasticModel>(args, kwargs, {});
        }))
      ;
}

}
