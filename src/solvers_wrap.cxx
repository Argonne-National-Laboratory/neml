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

  py::class_<SolverParameters, std::shared_ptr<SolverParameters>>(m,
                                                                  "SolverParameters")
      .def(py::init<double, double, int, bool, bool>())
      .def_readwrite("rtol", &SolverParameters::rtol)
      .def_readwrite("atol", &SolverParameters::atol)
      .def_readwrite("miter", &SolverParameters::miter)
      .def_readwrite("verbose", &SolverParameters::verbose)
      .def_readwrite("linesearch", &SolverParameters::linesearch)
      .def_readwrite("mline", &SolverParameters::mline)
    ;

  py::class_<Solvable, std::shared_ptr<Solvable>>(m, "Solvable")
      .def_property_readonly("nparams", &Solvable::nparams, "Number of variables in nonlinear equations.")
      .def("init_x",
           [](Solvable & m, TrialState & ts) -> py::array_t<double>
           {
            auto x = alloc_vec<double>(m.nparams());
            m.init_x(arr2ptr<double>(x), &ts);
            return x;
           }, "Initialize guess.")
      .def("RJ",
           [](Solvable & m, py::array_t<double, py::array::c_style> x, TrialState & ts) -> std::tuple<py::array_t<double>, py::array_t<double>>
           {
            auto R = alloc_vec<double>(m.nparams());
            auto J = alloc_mat<double>(m.nparams(), m.nparams());
            
            m.RJ(arr2ptr<double>(x), &ts, arr2ptr<double>(R), arr2ptr<double>(J));

            return std::make_tuple(R, J);
           }, "Residual and jacobian.")
      ;

  m.def("solve",
        [](std::shared_ptr<Solvable> system, TrialState & ts, double rtol,
           double atol, int miter, bool verbose, bool linesearch) -> py::array_t<double>
        {
          auto x = alloc_vec<double>(system->nparams());
          
          solve(system.get(), arr2ptr<double>(x), &ts, 
                          {rtol, atol, miter, verbose, linesearch});
          return x;
        }, "Solve a nonlinear system", 
        py::arg("solvable"), py::arg("trial_state"), py::arg("rtol") = 1.0e-6, py::arg("atol") = 1.0e-8,
        py::arg("miter") = 50,
        py::arg("verbose") = false, py::arg("linesearch") = false);

  py::class_<TestPower, Solvable, std::shared_ptr<TestPower>>(m, "TestPower")
      .def(py::init<double, double, double, double>())
      ;
}

} // namespace neml
