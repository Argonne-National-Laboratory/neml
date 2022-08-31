#include "pyhelp.h" // include first to avoid annoying redef warning

#include "walker.h"

#include "objects.h"
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
             m.kappa(arr2ptr<double>(eps), T, res);
             return res;
           }, "Rate sliding function")
      .def("dkappa",
           [](WalkerKremplSwitchRule & m, py::array_t<double,
              py::array::c_style> eps, double T) -> py::array_t<double>
           {
            auto f = alloc_vec<double>(6);
            m.dkappa(arr2ptr<double>(eps), T, arr2ptr<double>(f));
            return f;
           }, "Derivative of the kappa function wrt strain rate.")
      ;

  py::class_<SofteningModel, NEMLObject, std::shared_ptr<SofteningModel>>(m,
                                                              "SofteningModel")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<SofteningModel>(args, kwargs, {});
      }))
    .def("phi", &SofteningModel::phi)
    .def("dphi", &SofteningModel::dphi)
    ;

  py::class_<WalkerSofteningModel, SofteningModel, std::shared_ptr<WalkerSofteningModel>>(m,
                                                              "WalkerSofteningModel")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<WalkerSofteningModel>(
            args, kwargs, {"phi_0", "phi_1"});
      }))
    ;

  py::class_<ThermalScaling, NEMLObject, std::shared_ptr<ThermalScaling>>(m,
                                                              "ThermalScaling")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<ThermalScaling>(args, kwargs, {});
      }))
    .def("value", &ThermalScaling::value)
    ;

  py::class_<ArrheniusThermalScaling, ThermalScaling, std::shared_ptr<ArrheniusThermalScaling>>(m,
                                                              "ArrheniusThermalScaling")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<ArrheniusThermalScaling>(
            args, kwargs, {"Q", "R", "T_ref"});
      }))
    ;

  py::class_<ScalarInternalVariable::VariableState,
      std::shared_ptr<ScalarInternalVariable::VariableState>>(m,
                                                              "ScalarInternalVariableState")
    .def(py::init([](double h, double a, double adot, double D, Symmetric s, Symmetric g,
                     double T)
                  {
                    ScalarInternalVariable::VariableState state;
                    state.h = h;
                    state.a = a;
                    state.adot = adot;
                    state.D = D;
                    state.s = s;
                    state.g = g;
                    state.T = T;

                    return state;
                  }))
    .def_readwrite("h", &ScalarInternalVariable::VariableState::h)
    .def_readwrite("a", &ScalarInternalVariable::VariableState::a)
    .def_readwrite("adot", &ScalarInternalVariable::VariableState::adot)
    .def_readwrite("D", &ScalarInternalVariable::VariableState::D)
    .def_readwrite("s", &ScalarInternalVariable::VariableState::s)
    .def_readwrite("g", &ScalarInternalVariable::VariableState::g)
    .def_readwrite("T", &ScalarInternalVariable::VariableState::T)
    ;

  py::class_<ScalarInternalVariable, NEMLObject,
      std::shared_ptr<ScalarInternalVariable>>(m, "ScalarInternalVariable")
    .def("initial_value", &ScalarInternalVariable::initial_value)

    .def("ratep", &ScalarInternalVariable::ratep)
    .def("d_ratep_d_h", &ScalarInternalVariable::d_ratep_d_h)
    .def("d_ratep_d_a", &ScalarInternalVariable::d_ratep_d_a)
    .def("d_ratep_d_adot", &ScalarInternalVariable::d_ratep_d_adot)
    .def("d_ratep_d_D", &ScalarInternalVariable::d_ratep_d_D)
    .def("d_ratep_d_s", &ScalarInternalVariable::d_ratep_d_s)
    .def("d_ratep_d_g", &ScalarInternalVariable::d_ratep_d_g)

    .def("ratet", &ScalarInternalVariable::ratet)
    .def("d_ratet_d_h", &ScalarInternalVariable::d_ratet_d_h)
    .def("d_ratet_d_a", &ScalarInternalVariable::d_ratet_d_a)
    .def("d_ratet_d_adot", &ScalarInternalVariable::d_ratet_d_adot)
    .def("d_ratet_d_D", &ScalarInternalVariable::d_ratet_d_D)
    .def("d_ratet_d_s", &ScalarInternalVariable::d_ratet_d_s)
    .def("d_ratet_d_g", &ScalarInternalVariable::d_ratet_d_g)

    .def("rateT", &ScalarInternalVariable::rateT)
    .def("d_rateT_d_h", &ScalarInternalVariable::d_rateT_d_h)
    .def("d_rateT_d_a", &ScalarInternalVariable::d_rateT_d_a)
    .def("d_rateT_d_adot", &ScalarInternalVariable::d_rateT_d_adot)
    .def("d_rateT_d_D", &ScalarInternalVariable::d_rateT_d_D)
    .def("d_rateT_d_s", &ScalarInternalVariable::d_rateT_d_s)
    .def("d_rateT_d_g", &ScalarInternalVariable::d_rateT_d_g)
    ;

  py::class_<IsotropicHardening, ScalarInternalVariable,
      std::shared_ptr<IsotropicHardening>>(m, "IsotropicHardening")
      .def("set_scaling", &IsotropicHardening::set_scaling)
    ;

  py::class_<ConstantIsotropicHardening, IsotropicHardening,
      std::shared_ptr<ConstantIsotropicHardening>>(m,
                                                   "ConstantIsotropicHardening")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<ConstantIsotropicHardening>(
            args, kwargs, {});
      }))
    ;

  py::class_<WalkerIsotropicHardening, IsotropicHardening,
      std::shared_ptr<WalkerIsotropicHardening>>(m,
                                                   "WalkerIsotropicHardening")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<WalkerIsotropicHardening>(
            args, kwargs, {"r0", "Rinf", "R0", "r1", "r2"});
      }))
    ;

  py::class_<DragStress, ScalarInternalVariable,
      std::shared_ptr<DragStress>>(m, "DragStress")
      .def("set_scaling", &DragStress::set_scaling)
      .def("D_xi", &DragStress::D_xi)
      .def("D_0", &DragStress::D_0)
    ;

  py::class_<ConstantDragStress, DragStress,
      std::shared_ptr<ConstantDragStress>>(m, "ConstantDragStress")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<ConstantDragStress>(
            args, kwargs, {"value"});
      }))
    ;

  py::class_<WalkerDragStress, DragStress,
      std::shared_ptr<WalkerDragStress>>(m, "WalkerDragStress")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<WalkerDragStress>(
            args, kwargs, {"d0", "d1", "d2", "D_xi",
            "D_0", "softening"});
      }))
    ;

  py::class_<SymmetricInternalVariable::VariableState,
      std::shared_ptr<SymmetricInternalVariable::VariableState>>(m,
                                                              "SymmetricInternalVariableState")
    .def(py::init([](Symmetric h, double a, double adot, double D, Symmetric s, Symmetric g,
                     double T)
                  {
                    SymmetricInternalVariable::VariableState state;
                    state.h = h;
                    state.a = a;
                    state.adot = adot;
                    state.D = D;
                    state.s = s;
                    state.g = g;
                    state.T = T;

                    return state;
                  }))
    .def_readwrite("h", &SymmetricInternalVariable::VariableState::h)
    .def_readwrite("a", &SymmetricInternalVariable::VariableState::a)
    .def_readwrite("adot", &SymmetricInternalVariable::VariableState::adot)
    .def_readwrite("D", &SymmetricInternalVariable::VariableState::D)
    .def_readwrite("s", &SymmetricInternalVariable::VariableState::s)
    .def_readwrite("g", &SymmetricInternalVariable::VariableState::g)
    .def_readwrite("T", &SymmetricInternalVariable::VariableState::T)
    ;

  py::class_<SymmetricInternalVariable, NEMLObject,
      std::shared_ptr<SymmetricInternalVariable>>(m, "SymmetricInternalVariable")
    .def("initial_value", &SymmetricInternalVariable::initial_value)

    .def("ratep", &SymmetricInternalVariable::ratep)
    .def("d_ratep_d_h", &SymmetricInternalVariable::d_ratep_d_h)
    .def("d_ratep_d_a", &SymmetricInternalVariable::d_ratep_d_a)
    .def("d_ratep_d_adot", &SymmetricInternalVariable::d_ratep_d_adot)
    .def("d_ratep_d_D", &SymmetricInternalVariable::d_ratep_d_D)
    .def("d_ratep_d_s", &SymmetricInternalVariable::d_ratep_d_s)
    .def("d_ratep_d_g", &SymmetricInternalVariable::d_ratep_d_g)

    .def("ratet", &SymmetricInternalVariable::ratet)
    .def("d_ratet_d_h", &SymmetricInternalVariable::d_ratet_d_h)
    .def("d_ratet_d_a", &SymmetricInternalVariable::d_ratet_d_a)
    .def("d_ratet_d_adot", &SymmetricInternalVariable::d_ratet_d_adot)
    .def("d_ratet_d_D", &SymmetricInternalVariable::d_ratet_d_D)
    .def("d_ratet_d_s", &SymmetricInternalVariable::d_ratet_d_s)
    .def("d_ratet_d_g", &SymmetricInternalVariable::d_ratet_d_g)

    .def("rateT", &SymmetricInternalVariable::rateT)
    .def("d_rateT_d_h", &SymmetricInternalVariable::d_rateT_d_h)
    .def("d_rateT_d_a", &SymmetricInternalVariable::d_rateT_d_a)
    .def("d_rateT_d_adot", &SymmetricInternalVariable::d_rateT_d_adot)
    .def("d_rateT_d_D", &SymmetricInternalVariable::d_rateT_d_D)
    .def("d_rateT_d_s", &SymmetricInternalVariable::d_rateT_d_s)
    .def("d_rateT_d_g", &SymmetricInternalVariable::d_rateT_d_g)
    ;

  py::class_<KinematicHardening, SymmetricInternalVariable,
      std::shared_ptr<KinematicHardening>>(m, "KinematicHardening")
      .def("set_scaling", &KinematicHardening::set_scaling)
    ;

  py::class_<FAKinematicHardening, KinematicHardening,
      std::shared_ptr<FAKinematicHardening>>(m, "FAKinematicHardening")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<FAKinematicHardening>(
            args, kwargs, {"c", "g"});
      }))
    ;

  py::class_<WalkerKinematicHardening, KinematicHardening,
      std::shared_ptr<WalkerKinematicHardening>>(m, "WalkerKinematicHardening")
    .def(py::init([](py::args args, py::kwargs kwargs)
      {
        return create_object_python<WalkerKinematicHardening>(
            args, kwargs, {"c0", "c1", "c2", "l0", "l1", "l",
            "b0", "x0", "x1", "softening"});
      }))
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
          History h = m.gather_blank_history_();
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
          History res = m.gather_blank_history_().derivative<Symmetric>();
          m.dg_da(state, res);
          return res;
         }, "The derivative of the flow direction wrt history")

    .def("h_wrap",
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_();
          m.h(state, res);
          return res;
         }, "The history proportional to the scalar inelastic strain rate")
    .def("dh_ds_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_().derivative<Symmetric>();
          m.dh_ds(state, res);
          return res;
         }, "The derivative of the history ptt inelastic strain rate wrt stress")
    .def("dh_da_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_().history_derivative(m.gather_blank_history_());
          m.dh_da(state, res);
          return res;
         }, "The derivative of the history ptt inelastic strain rate wrt history")

    .def("h_time_wrap",
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_();
          m.h_time(state, res);
          return res;
         }, "The history proportional to the time rate")
    .def("dh_ds_time_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_().derivative<Symmetric>();
          m.dh_ds_time(state, res);
          return res;
         }, "The derivative of the history ptt time rate wrt stress")
    .def("dh_da_time_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_().history_derivative(m.gather_blank_history_());
          m.dh_da_time(state, res);
          return res;
         }, "The derivative of the history ptt time rate wrt history")

    .def("h_temp_wrap",
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_();
          m.h_temp(state, res);
          return res;
         }, "The history proportional to the temperature rate")
    .def("dh_ds_temp_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_().derivative<Symmetric>();
          m.dh_ds_temp(state, res);
          return res;
         }, "The derivative of the history ptt temperature rate wrt stress")
    .def("dh_da_temp_wrap", 
         [](WrappedViscoPlasticFlowRule & m, State & state) -> History
         {
          History res = m.gather_blank_history_().history_derivative(m.gather_blank_history_());
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

  py::class_<WalkerFlowRule, WrappedViscoPlasticFlowRule,
      std::shared_ptr<WalkerFlowRule>>(m, "WalkerFlowRule")
    .def(py::init([](py::args args, py::kwargs kwargs)
          {
            return create_object_python<WalkerFlowRule>(args, kwargs,
                    {"eps0", "softening", "scaling", "n", "k", "m", "R",
                    "D", "X"});
          }))
    ;
}

} // namespace neml
