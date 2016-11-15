#include "solvers.h"

#include "pyhelp.h"
#include "nemlerror.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

namespace neml {

PYBIND11_PLUGIN(solvers) {
  py::module m("solvers", "Nonlinear solvers.");

  py::class_<Solvable, std::shared_ptr<Solvable>>(m, "Solvable")
      .def_property_readonly("nparams", &Solvable::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](Solvable & m) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x));
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](Solvable & m, py::array_t<double, py::array::c_style> x) -> std::tuple<py::array_t<double>, py::array_t<double>>
           {
            auto R = alloc_vec<double>(m.nparams());
            auto J = alloc_mat<double>(m.nparams(), m.nparams());
            
            int ier = m.RJ(arr2ptr<double>(x), arr2ptr<double>(R), arr2ptr<double>(J));
            py_error(ier);

            return std::make_tuple(R, J);
           }, "Residual and jacobian.")
      ;

  py::class_<TestRosenbrock, std::shared_ptr<TestRosenbrock>>(m, "TestRosenbrock", py::base<Solvable>())
      .def(py::init<size_t>(), py::arg("N"))
      ;

  m.def("solve",
        [](std::shared_ptr<Solvable> system, double tol, int miter, bool verbose) -> py::array_t<double>
        {
          auto x = alloc_vec<double>(system->nparams());
          
          int ier = solve(system, arr2ptr<double>(x), tol, miter, verbose);
          py_error(ier);

          return x;
        }, "Solve a nonlinear system", 
        py::arg("solvable"), py::arg("tol") = 1.0e-8,
        py::arg("miter") = 50,
        py::arg("verbose") = false);

  return m.ptr();

}

} // namespace neml
