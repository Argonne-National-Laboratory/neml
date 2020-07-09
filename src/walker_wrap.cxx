#include "pyhelp.h" // include first to avoid annoying redef warning

#include "walker.h"

#include "nemlerror.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {
PYBIND11_MODULE(walker, m) {
  py::module::import("neml.objects");
  py::module::import("neml.general_flow");
  py::module::import("neml.visco_flow");

  m.doc() = "Objects specific to Walker's A617 model";

  py::class_<WalkerKremplSwitchRule, GeneralFlowRule,
      std::shared_ptr<WalkerKremplSwitchRule>>(m, "WalkerKremplSwitchRule")
      .def(py::init([](py::args args, py::kwargs kwargs)
        {
          return create_object_python<WalkerKremplSwitchRule>(args, kwargs,
                                                              {"elastic",
                                                              "flow", "lambda",
                                                              "eps_ref"});
        }))
      .def("kappa",
           [](WalkerKremplSwitchRule & m, py::array_t<double,
              py::array::c_style> eps, double T) -> double
           {
             double res;
             int ier = m.kappa(arr2ptr<double>(eps), T, res);
             py_error(ier);
             return res;
           }, "Rate sliding function")
      .def("dkappa",
           [](WalkerKremplSwitchRule & m, py::array_t<double,
              py::array::c_style> eps, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            int ier = m.dkappa(arr2ptr<double>(eps), T, arr2ptr<double>(f));
            py_error(ier);
            return f;
           }, "Derivative of the kappa function wrt strain rate.")
      ;

  py::class_<State, std::shared_ptr<State>>(m, "State")
    .def(py::init<Symmetric, History, double>())
    .def_readwrite("S", &State::S)
    .def_readwrite("h", &State::h)
    .def_readwrite("T", &State::T)
    ;

  py::class_<WrappedViscoPlasticFlowRule, ViscoPlasticFlowRule,
      std::shared_ptr<WrappedViscoPlasticFlowRule>>(m,
                                                    "WrappedViscoPlasticFlowRule")
    .def("populate_hist",
         [](WrappedViscoPlasticFlowRule & m) -> History
         {
          History h;
          m.populate_hist(h);
          return h;
         }, "Return a blank history")
    .def("initialize_hist", 
         [](WrappedViscoPlasticFlowRule & m) -> History
         {
          History h;
          m.populate_hist(h);
          m.initialize_hist(h);
          return h;
         }, "Initialize the history, including setting the initial conditions.")
    .def("y_wrap",
         [](WrappedViscoPlasticFlowRule & m, State & state) -> double
         {
          double yv;
          m.y(state, yv);
          return yv;
         }, "The scalar inelastic strain rate")
    .def("dy_ds_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> Symmetric
         {
          Symmetric s;
          m.dy_ds(state, s);
          return s;
         }, "The derivative of the scalar inelastic strain rate wrt stress")
    .def("dy_da_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History h = m.create_blank_hist_();
          m.dy_da(state, h);
          return h;
         }, "The derivative of the scalar inelastic strain rate wrt history")

    .def("g_wrap",
         [](WrappedViscoPlasticFlowRule & m, State & state) -> Symmetric
         {
          Symmetric res;
          m.g(state, res);
          return res;
         }, "The flow direction")
    .def("dg_ds_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> SymSymR4
         {
          SymSymR4 res;
          m.dg_ds(state, res);
          return res;
         }, "The derivative of the flow direction wrt stress")
    .def("dg_da_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_().derivative<Symmetric>();
          m.dg_da(state, res);
          return res;
         }, "The derivative of the flow direction wrt history")

    .def("h_wrap",
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_();
          m.h(state, res);
          return res;
         }, "The history proportional to the scalar inelastic strain rate")
    .def("dh_ds_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_().derivative<Symmetric>();
          m.dh_ds(state, res);
          return res;
         }, "The derivative of the history ptt inelastic strain rate wrt stress")
    .def("dh_da_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_().history_derivative(m.create_blank_hist_());
          m.dh_da(state, res);
          return res;
         }, "The derivative of the history ptt inelastic strain rate wrt history")

    .def("h_time_wrap",
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_();
          m.h_time(state, res);
          return res;
         }, "The history proportional to the time rate")
    .def("dh_ds_time_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_().derivative<Symmetric>();
          m.dh_ds_time(state, res);
          return res;
         }, "The derivative of the history ptt time rate wrt stress")
    .def("dh_da_time_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_().history_derivative(m.create_blank_hist_());
          m.dh_da_time(state, res);
          return res;
         }, "The derivative of the history ptt time rate wrt history")

    .def("h_temp_wrap",
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_();
          m.h_temp(state, res);
          return res;
         }, "The history proportional to the temperature rate")
    .def("dh_ds_temp_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_().derivative<Symmetric>();
          m.dh_ds_temp(state, res);
          return res;
         }, "The derivative of the history ptt temperature rate wrt stress")
    .def("dh_da_temp_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.create_blank_hist_().history_derivative(m.create_blank_hist_());
          m.dh_da_temp(state, res);
          return res;
         }, "The derivative of the history ptt temperature rate wrt history")
    ;

  py::class_<TestFlowRule, WrappedViscoPlasticFlowRule,
      std::shared_ptr<TestFlowRule>>(m, "TestFlowRule")
    .def(py::init([](py::args args, py::kwargs kwargs)
                  {
                    return create_object_python<TestFlowRule>(args, kwargs,
                                                              {"eps0", "D", "n",
                                                              "s0", "K"});
                  }))
    ;
}

} // namespace neml
