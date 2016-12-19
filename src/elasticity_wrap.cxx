#include "elasticity.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(elasticity) {
  py::module m("elasticity", "Models of elasticity.");

  py::class_<ShearModulus, std::shared_ptr<ShearModulus>>(m, "ShearModulus")
      .def("modulus", &ShearModulus::modulus, "Modulus as a function of temperature.")
      ;

  py::class_<ConstantShearModulus, std::shared_ptr<ConstantShearModulus>>(m, "ConstantShearModulus", py::base<ShearModulus>())
      .def(py::init<double>(), py::arg("mu"))
      .def_property_readonly("mu", &ConstantShearModulus::mu, "Constant modulus value")
      ;

  py::class_<PolyShearModulus, std::shared_ptr<PolyShearModulus>>(m, "PolyShearModulus", py::base<ShearModulus>())
      .def(py::init<std::vector<double>>(), py::arg("coefs"))
      .def_property_readonly("n", &PolyShearModulus::n, "Number of coefs")
      .def_property_readonly("coefs",
                             [](const PolyShearModulus& m) -> py::array_t<double>
                             {
                              auto cv = alloc_vec<double>(m.n());
                              std::copy(m.coefs().begin(), m.coefs().end(), arr2ptr<double>(cv));
                              return cv;
                             }, "Polynomial coefficients.")
      ;

  py::class_<BulkModulus, std::shared_ptr<BulkModulus>>(m, "BulkModulus")
      .def("modulus", &BulkModulus::modulus, "Modulus as a function of temperature.")
      ;

  py::class_<ConstantBulkModulus, std::shared_ptr<ConstantBulkModulus>>(m, "ConstantBulkModulus", py::base<BulkModulus>())
      .def(py::init<double>(), py::arg("K"))
      .def_property_readonly("K", &ConstantBulkModulus::K, "Constant modulus value")
      ;

  py::class_<PolyBulkModulus, std::shared_ptr<PolyBulkModulus>>(m, "PolyBulkModulus", py::base<BulkModulus>())
      .def(py::init<std::vector<double>>(), py::arg("coefs"))
      .def_property_readonly("n", &PolyBulkModulus::n, "Number of coefs")
      .def_property_readonly("coefs",
                             [](const PolyBulkModulus& m) -> py::array_t<double>
                             {
                              auto cv = alloc_vec<double>(m.n());
                              std::copy(m.coefs().begin(), m.coefs().end(), arr2ptr<double>(cv));
                              return cv;
                             }, "Polynomial coefficients.")
      ;

  py::class_<LinearElasticModel, std::shared_ptr<LinearElasticModel>>(m, "LinearElasticModel")
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
      ;

  py::class_<IsotropicLinearElasticModel, std::shared_ptr<IsotropicLinearElasticModel>>(m, "IsotropicLinearElasticModel", py::base<LinearElasticModel>())
      .def(py::init<std::shared_ptr<ShearModulus>, std::shared_ptr<BulkModulus>>(),
           py::arg("shear"), py::arg("bulk"))
      ;

  return m.ptr();
}


}
