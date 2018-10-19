#include "pyhelp.h" // include first to avoid annoying redef warning

#include "solvers.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(solvers, m) {
  m.doc() = "Nonlinear solvers and wrappers to nonlinear solver libraries.";

  py::class_<TrialState>(m, "TrialState")
      .def(py::init<>())
      ;

  py::class_<Solvable, std::shared_ptr<Solvable>>(m, "Solvable")
      .def_property_readonly("nparams", &Solvable::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](Solvable & m, TrialState & ts) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            int ier = m.init_x(arr2ptr<double>(x), &ts);
            py_error(ier);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](Solvable & m, py::array_t<double, py::array::c_style> x, TrialState & ts) -> std::tuple<py::array_t<double>, py::array_t<double>>
           {
            auto R = alloc_vec<double>(m.nparams());
            auto J = alloc_mat<double>(m.nparams(), m.nparams());
            
            int ier = m.RJ(arr2ptr<double>(x), &ts, arr2ptr<double>(R), arr2ptr<double>(J));
            py_error(ier);

            return std::make_tuple(R, J);
           }, "Residual and jacobian.")
      ;

  m.def("solve",
        [](std::shared_ptr<Solvable> system, TrialState & ts, double tol, int miter, bool verbose) -> py::array_t<double>
        {
          auto x = alloc_vec<double>(system->nparams());
          
          int ier = solve(system.get(), arr2ptr<double>(x), &ts, tol, miter, verbose);
          py_error(ier);

          return x;
        }, "Solve a nonlinear system", 
        py::arg("solvable"), py::arg("trial_state"), py::arg("tol") = 1.0e-8,
        py::arg("miter") = 50,
        py::arg("verbose") = false);
}

} // namespace neml
