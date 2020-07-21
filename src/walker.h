#pragma once

#include "history.h"
#include "math/tensors.h"

#include "general_flow.h"
#include "visco_flow.h"

#include "internalvariable.h" // Mostly ABCs for the "fancy" history classes

#include "windows.h"

namespace neml {

/// WalkerKremplSwitchRule
//    Switches between rate independent and rate dependent flow based on 
//    a parameter "gamma"
//    This implementation ignores any time or temperature rate dependent
//    *flow* terms in the flow rule (it includes those terms in the hardening)
class NEML_EXPORT WalkerKremplSwitchRule : public GeneralFlowRule {
 public:
  /// Parameters: elastic model and a viscoplastic flow rule
  WalkerKremplSwitchRule(std::shared_ptr<LinearElasticModel> elastic,
                         std::shared_ptr<ViscoPlasticFlowRule> flow,
                         std::shared_ptr<Interpolate> lambda,
                         double eps0);

  /// String type for the object system
  static std::string type();
  /// Return default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from parameter set
  static ParameterSet parameters();

  /// Number of history variables
  virtual size_t nhist() const;
  /// Initialize history
  virtual int init_hist(double * const h);

  /// Stress rate
  virtual int s(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const sdot);
  /// Partial of stress rate wrt stress
  virtual int ds_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);
  /// Partial of stress rate wrt history
  virtual int ds_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);
  /// Partial of stress rate wrt strain rate
  virtual int ds_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_sdot);

  /// History rate
  virtual int a(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const adot);
  /// Partial of history rate wrt stress
  virtual int da_ds(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);
  /// Partial of history rate wrt history
  virtual int da_da(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);
  /// Partial of history rate wrt strain rate
  virtual int da_de(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double * const d_adot);

  /// The implementation needs to define inelastic dissipation
  virtual int work_rate(const double * const s, const double * const alpha,
                const double * const edot, double T,
                double Tdot,
                double & p_rate);

  /// The implementation needs to define elastic strain
  virtual int elastic_strains(const double * const s_np1, double T_np1,
                              double * const e_np1) const;

  /// Set a new elastic model
  virtual int set_elastic_model(std::shared_ptr<LinearElasticModel> emodel);
  
  /// The kappa function controlling rate sensitivity (public for testing)
  int kappa(const double * const edot, double T, double & kap);

  /// Derivative of kappa wrt the strain rate (public for testing)
  int dkappa(const double * const edot, double T, double * const dkap);

 private:
  std::shared_ptr<LinearElasticModel> elastic_;
  std::shared_ptr<ViscoPlasticFlowRule> flow_;
  std::shared_ptr<Interpolate> lambda_;
  const double eps0_;
};

static Register<WalkerKremplSwitchRule> regWalkerKremplSwitchRule;

/// Softening/tertiary creep models
class NEML_EXPORT SofteningModel: public NEMLObject {
 public:
  SofteningModel();

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// Softening function
  virtual double phi(double alpha, double T) const;
  /// Derivative of softening function wrt alpha
  virtual double dphi(double alpha, double T) const;
};

static Register<SofteningModel> regSoftening;

/// The simple softening model Walker actually uses
class NEML_EXPORT WalkerSofteningModel: public SofteningModel {
 public:
  WalkerSofteningModel(std::shared_ptr<Interpolate> phi0,
                       std::shared_ptr<Interpolate> phi1);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// Softening function
  virtual double phi(double alpha, double T) const;
  /// Derivative of softening function wrt alpha
  virtual double dphi(double alpha, double T) const;

 private:
  std::shared_ptr<Interpolate> phi_0_;
  std::shared_ptr<Interpolate> phi_1_;
};

static Register<WalkerSofteningModel> regWalkerSoftening;

/// Thermal rate scaling models
class ThermalScaling: public NEMLObject {
 public:
  ThermalScaling();

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// The thermal ratio
  virtual double value(double T) const;
};

static Register<ThermalScaling> regThermalScaling;

/// Walker's actual thermal scaling model
class ArrheniusThermalScaling: public ThermalScaling {
 public:
  ArrheniusThermalScaling(std::shared_ptr<Interpolate> Q,
                          double R, double Tref);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  /// The thermal ratio
  virtual double value(double T) const;

 private:
  /// The actual Arrhenius function
  double arr_(double T) const;

 private:
  std::shared_ptr<Interpolate> Q_;
  double R_;
  double T_ref_;
};

static Register<ArrheniusThermalScaling> regArrheniusThermalScaling;

class IsotropicHardening: public ScalarInternalVariable {
 public:
  IsotropicHardening(std::string name, 
                     std::shared_ptr<ThermalScaling> scale);

  void set_scaling(std::shared_ptr<ThermalScaling> scale) {scale_ = scale;};

  virtual double ratet(VariableState & state);
  virtual double d_ratet_d_h(VariableState & state);
  virtual double d_ratet_d_a(VariableState & state);
  virtual double d_ratet_d_adot(VariableState & state);
  virtual double d_ratet_d_D(VariableState & state);
  virtual Symmetric d_ratet_d_s(VariableState & state);
  virtual Symmetric d_ratet_d_g(VariableState & state);

  virtual double rateT(VariableState & state);
  virtual double d_rateT_d_h(VariableState & state);
  virtual double d_rateT_d_a(VariableState & state);
  virtual double d_rateT_d_adot(VariableState & state);
  virtual double d_rateT_d_D(VariableState & state);
  virtual Symmetric d_rateT_d_s(VariableState & state);
  virtual Symmetric d_rateT_d_g(VariableState & state);

 protected:
  std::shared_ptr<ThermalScaling> scale_;
};

class ConstantIsotropicHardening: public IsotropicHardening {
 public:
  ConstantIsotropicHardening(std::shared_ptr<ThermalScaling> scale = 
                             std::make_shared<ThermalScaling>());

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  virtual double initial_value();

  virtual double ratep(VariableState & state);
  virtual double d_ratep_d_h(VariableState & state);
  virtual double d_ratep_d_a(VariableState & state);
  virtual double d_ratep_d_adot(VariableState & state);
  virtual double d_ratep_d_D(VariableState & state);
  virtual Symmetric d_ratep_d_s(VariableState & state);
  virtual Symmetric d_ratep_d_g(VariableState & state);
};

static Register<ConstantIsotropicHardening> regConstantIsotropicHardening;

/// The actual hardening model used in Walker's A617 viscoplastic model
class WalkerIsotropicHardening: public IsotropicHardening {
 public:
  WalkerIsotropicHardening(std::shared_ptr<Interpolate> r0,
                           std::shared_ptr<Interpolate> Rinf,
                           std::shared_ptr<Interpolate> R0,
                           std::shared_ptr<Interpolate> r1,
                           std::shared_ptr<Interpolate> r2,
                           std::shared_ptr<ThermalScaling> scale = 
                           std::make_shared<ThermalScaling>());

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  virtual double initial_value();

  virtual double ratep(VariableState & state);
  virtual double d_ratep_d_h(VariableState & state);
  virtual double d_ratep_d_a(VariableState & state);
  virtual double d_ratep_d_adot(VariableState & state);
  virtual double d_ratep_d_D(VariableState & state);
  virtual Symmetric d_ratep_d_s(VariableState & state);
  virtual Symmetric d_ratep_d_g(VariableState & state);

  virtual double ratet(VariableState & state);
  virtual double d_ratet_d_h(VariableState & state);
  virtual double d_ratet_d_a(VariableState & state);
  virtual double d_ratet_d_adot(VariableState & state);
  virtual double d_ratet_d_D(VariableState & state);
  virtual Symmetric d_ratet_d_s(VariableState & state);
  virtual Symmetric d_ratet_d_g(VariableState & state);

 private:
  std::shared_ptr<Interpolate> r0_;
  std::shared_ptr<Interpolate> Rinf_;
  std::shared_ptr<Interpolate> R0_;
  std::shared_ptr<Interpolate> r1_;
  std::shared_ptr<Interpolate> r2_;
};

static Register<WalkerIsotropicHardening> regWalkerIsotropicHardening;

class DragStress: public ScalarInternalVariable {
 public:
  DragStress(std::string name, 
             std::shared_ptr<ThermalScaling> scale);

  void set_scaling(std::shared_ptr<ThermalScaling> scale) {scale_ = scale;};
  
  /// Report the value of D_xi
  virtual double D_xi(double T) = 0;
  /// Report the value of D_0
  virtual double D_0(double T) = 0;
  
  /// Makes no sense in this context!
  virtual double d_ratep_d_D(VariableState & state);

  virtual double ratet(VariableState & state);
  virtual double d_ratet_d_h(VariableState & state);
  virtual double d_ratet_d_a(VariableState & state);
  virtual double d_ratet_d_adot(VariableState & state);
  virtual double d_ratet_d_D(VariableState & state);
  virtual Symmetric d_ratet_d_s(VariableState & state);
  virtual Symmetric d_ratet_d_g(VariableState & state);

  virtual double rateT(VariableState & state);
  virtual double d_rateT_d_h(VariableState & state);
  virtual double d_rateT_d_a(VariableState & state);
  virtual double d_rateT_d_adot(VariableState & state);
  virtual double d_rateT_d_D(VariableState & state);
  virtual Symmetric d_rateT_d_s(VariableState & state);
  virtual Symmetric d_rateT_d_g(VariableState & state);

 protected:
  std::shared_ptr<ThermalScaling> scale_;
};

class ConstantDragStress: public DragStress {
 public:
  ConstantDragStress(double value,
                     std::shared_ptr<ThermalScaling> scale = 
                     std::make_shared<ThermalScaling>());

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  virtual double initial_value();

  /// Report the value of D_xi
  virtual double D_xi(double T);
  /// Report the value of D_0
  virtual double D_0(double T);

  virtual double ratep(VariableState & state);
  virtual double d_ratep_d_h(VariableState & state);
  virtual double d_ratep_d_a(VariableState & state);
  virtual double d_ratep_d_adot(VariableState & state);
  virtual Symmetric d_ratep_d_s(VariableState & state);
  virtual Symmetric d_ratep_d_g(VariableState & state);

 private:
  double value_;
};

static Register<ConstantDragStress> regConstantDragStress;

class WalkerDragStress: public DragStress {
 public:
  WalkerDragStress(std::shared_ptr<Interpolate> d0,
                   std::shared_ptr<Interpolate> d1,
                   std::shared_ptr<Interpolate> d2,
                   std::shared_ptr<Interpolate> D_xi,
                   double D_0,
                   std::shared_ptr<SofteningModel> softening,
                   std::shared_ptr<ThermalScaling> scale = 
                   std::make_shared<ThermalScaling>());

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();
  
  /// Initial value of drag stress
  virtual double initial_value();

  /// Report the value of D_xi
  virtual double D_xi(double T);
  /// Report the value of D_0
  virtual double D_0(double T);

  virtual double ratep(VariableState & state);
  virtual double d_ratep_d_h(VariableState & state);
  virtual double d_ratep_d_a(VariableState & state);
  virtual double d_ratep_d_adot(VariableState & state);
  virtual Symmetric d_ratep_d_s(VariableState & state);
  virtual Symmetric d_ratep_d_g(VariableState & state);

  virtual double ratet(VariableState & state);
  virtual double d_ratet_d_h(VariableState & state);
  virtual double d_ratet_d_a(VariableState & state);
  virtual double d_ratet_d_adot(VariableState & state);
  virtual Symmetric d_ratet_d_s(VariableState & state);
  virtual Symmetric d_ratet_d_g(VariableState & state);

 private:
  std::shared_ptr<Interpolate> d0_, d1_, d2_, D_xi_;
  double D_0_;
  std::shared_ptr<SofteningModel> softening_;
};

static Register<WalkerDragStress> regWalkerDragStress;

class KinematicHardening: public SymmetricInternalVariable {
 public:
  KinematicHardening(std::string name, 
                     std::shared_ptr<ThermalScaling> scale);

  void set_scaling(std::shared_ptr<ThermalScaling> scale) {scale_ = scale;};

  virtual Symmetric ratet(VariableState & state);
  virtual SymSymR4 d_ratet_d_h(VariableState & state);
  virtual Symmetric d_ratet_d_a(VariableState & state);
  virtual Symmetric d_ratet_d_adot(VariableState & state);
  virtual Symmetric d_ratet_d_D(VariableState & state);
  virtual SymSymR4 d_ratet_d_s(VariableState & state);
  virtual SymSymR4 d_ratet_d_g(VariableState & state);

  virtual Symmetric rateT(VariableState & state);
  virtual SymSymR4 d_rateT_d_h(VariableState & state);
  virtual Symmetric d_rateT_d_a(VariableState & state);
  virtual Symmetric d_rateT_d_adot(VariableState & state);
  virtual Symmetric d_rateT_d_D(VariableState & state);
  virtual SymSymR4 d_rateT_d_s(VariableState & state);
  virtual SymSymR4 d_rateT_d_g(VariableState & state);

 protected:
  std::shared_ptr<ThermalScaling> scale_;
};

/// Standard Frederick-Armstrong hardening
class FAKinematicHardening: public KinematicHardening {
 public:
  FAKinematicHardening(std::shared_ptr<Interpolate> c,
                       std::shared_ptr<Interpolate> g,
                       std::shared_ptr<ThermalScaling> scale);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  virtual Symmetric initial_value();

  virtual Symmetric ratep(VariableState & state);
  virtual SymSymR4 d_ratep_d_h(VariableState & state);
  virtual Symmetric d_ratep_d_a(VariableState & state);
  virtual Symmetric d_ratep_d_adot(VariableState & state);
  virtual Symmetric d_ratep_d_D(VariableState & state);
  virtual SymSymR4 d_ratep_d_s(VariableState & state);
  virtual SymSymR4 d_ratep_d_g(VariableState & state);

 private:
  std::shared_ptr<Interpolate> c_, g_;
};

static Register<FAKinematicHardening> regFAKinematicHardening;

/// Walker's kinematic hardening model
class WalkerKinematicHardening: public KinematicHardening {
 public:
  WalkerKinematicHardening(std::shared_ptr<Interpolate> c0,
                           std::shared_ptr<Interpolate> c1,
                           std::shared_ptr<Interpolate> c2,
                           std::shared_ptr<Interpolate> l0,
                           std::shared_ptr<Interpolate> l1,
                           std::shared_ptr<Interpolate> l,
                           std::shared_ptr<Interpolate> b0,
                           std::shared_ptr<Interpolate> x0,
                           std::shared_ptr<Interpolate> x1,
                           std::shared_ptr<SofteningModel> softening,
                           std::shared_ptr<ThermalScaling> scale);

  /// String type for the object system
  static std::string type();
  /// Initialize from a parameter set
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Return default parameters
  static ParameterSet parameters();

  virtual Symmetric initial_value();

  virtual Symmetric ratep(VariableState & state);
  virtual SymSymR4 d_ratep_d_h(VariableState & state);
  virtual Symmetric d_ratep_d_a(VariableState & state);
  virtual Symmetric d_ratep_d_adot(VariableState & state);
  virtual Symmetric d_ratep_d_D(VariableState & state);
  virtual SymSymR4 d_ratep_d_s(VariableState & state);
  virtual SymSymR4 d_ratep_d_g(VariableState & state);

  virtual Symmetric ratet(VariableState & state);
  virtual SymSymR4 d_ratet_d_h(VariableState & state);
  virtual Symmetric d_ratet_d_a(VariableState & state);
  virtual Symmetric d_ratet_d_adot(VariableState & state);
  virtual Symmetric d_ratet_d_D(VariableState & state);
  virtual SymSymR4 d_ratet_d_s(VariableState & state);
  virtual SymSymR4 d_ratet_d_g(VariableState & state);
 
 private:
  double c_(VariableState & state);
  double dc_(VariableState & state);
  double L_(VariableState & state);
  double dL_(VariableState & state);
  Symmetric n_(VariableState & state);
  Symmetric b_(VariableState & state);
  SymSymR4 dN_(VariableState & state);
  SymSymR4 db_ds_(VariableState & state);
  SymSymR4 db_dx_(VariableState & state);

 private:
  std::shared_ptr<Interpolate> c0_, c1_, c2_, l0_, l1_, l_, b0_, x0_, x1_;
  std::shared_ptr<SofteningModel> softening_;
};

static Register<WalkerKinematicHardening> regWalkerKinematicHardening;

/// Helper struct for the below
struct State {
  State(Symmetric S, History h, double T) :
      S(S), h(h), T(T) {};
  Symmetric S;
  History h;
  double T;
};

/// Wrapper between ViscoPlasticFlowRule and a version using the "fancy" objects
class NEML_EXPORT WrappedViscoPlasticFlowRule : public ViscoPlasticFlowRule {
 public:
  WrappedViscoPlasticFlowRule();

  /// Populate a history object
  virtual void populate_hist(History & h) const = 0;
  virtual void initialize_hist(History & h) const = 0;

  /// Number of history variables (from the hardening model)
  virtual size_t nhist() const;
  /// Initialize history at time zero
  virtual int init_hist(double * const h) const;

  // Scalar inelastic strain rate
  virtual int y(const double* const s, const double* const alpha, double T,
                double & yv) const;
  virtual void y(const State & state, double & res) const = 0;
  /// Derivative of y wrt stress
  virtual int dy_ds(const double* const s, const double* const alpha, double T,
                double * const dyv) const;
  virtual void dy_ds(const State & state, Symmetric & res) const = 0;
  /// Derivative of y wrt history
  virtual int dy_da(const double* const s, const double* const alpha, double T,
                double * const dyv) const;
  virtual void dy_da(const State & state, History & res) const = 0;

  /// Flow rule proportional to the scalar strain rate
  virtual int g(const double * const s, const double * const alpha, double T,
                double * const gv) const;
  virtual void g(const State & state, Symmetric & res) const = 0;
  /// Derivative of g wrt stress
  virtual int dg_ds(const double * const s, const double * const alpha, double T,
                double * const dgv) const;
  virtual void dg_ds(const State & state, SymSymR4 & res) const = 0;
  /// Derivative of g wrt history
  virtual int dg_da(const double * const s, const double * const alpha, double T,
               double * const dgv) const;
  virtual void dg_da(const State & state, History & res) const = 0;

  /// Hardening rule proportional to the scalar strain rate
  virtual int h(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual void h(const State & state, History & res) const = 0;
  /// Derivative of h wrt stress
  virtual int dh_ds(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual void dh_ds(const State & state, History & res) const = 0;
  /// Derivative of h wrt history
  virtual int dh_da(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual void dh_da(const State & state, History & res) const = 0;

  /// Hardening rule proportional to time
  virtual int h_time(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual void h_time(const State & state, History & res) const;
  /// Derivative of h_time wrt stress
  virtual int dh_ds_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual void dh_ds_time(const State & state, History & res) const;
  /// Derivative of h_time wrt history
  virtual int dh_da_time(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual void dh_da_time(const State & state, History & res) const;

  /// Hardening rule proportional to temperature rate
  virtual int h_temp(const double * const s, const double * const alpha, double T,
                double * const hv) const;
  virtual void h_temp(const State & state, History & res) const;
  /// Derivative of h_temp wrt stress
  virtual int dh_ds_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual void dh_ds_temp(const State & state, History & res) const;
  /// Derivative of h_temp wrt history
  virtual int dh_da_temp(const double * const s, const double * const alpha, double T,
                double * const dhv) const;
  virtual void dh_da_temp(const State & state, History & res) const;
  
  /// Blank history
  History blank_hist_() const;

  /// Special function for wrappers...
  History create_blank_hist_() const;
  
  /// Blank derivative of a history
  template <class T>
  History blank_derivative_() const
  {
    return blank_hist_().derivative<T>(); 
  }

 private:
  /// Make a state object
  State make_state_(const double * const s, const double * const alpha, double
                   T) const;

  /// Make a history object
  History gather_hist_(double * const h) const;
  History gather_hist_(const double * const h) const;

  /// Initialized derivative
  template <class T>
  History gather_derivative_(double * const h) const
  {
    History hv = blank_derivative_<T>();
    hv.set_data(h);
    return hv;
  }

 protected:
  History stored_hist_;
};

/// Initialized derivative wrt history
template <>
History WrappedViscoPlasticFlowRule::gather_derivative_<History>(double * const h) const
{
  History hv = blank_hist_().history_derivative(blank_hist_());
  hv.set_data(h);
  return hv;
}

/// Test implementation of a simple flow rule
class TestFlowRule: public WrappedViscoPlasticFlowRule
{
 public:
  TestFlowRule(double eps0, double D, double n, double s0, double K);

  /// String type for the object system
  static std::string type();
  /// Return default parameters
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);
  /// Initialize from parameter set
  static ParameterSet parameters();

  /// Populate a history object
  virtual void populate_hist(History & h) const;
  virtual void initialize_hist(History & h) const;

  // Scalar inelastic strain rate
  virtual void y(const State & state, double & res) const;
  /// Derivative of y wrt stress
  virtual void dy_ds(const State & state, Symmetric & res) const;
  /// Derivative of y wrt history
  virtual void dy_da(const State & state, History & res) const;

  /// Flow rule proportional to the scalar strain rate
  virtual void g(const State & state, Symmetric & res) const;
  /// Derivative of g wrt stress
  virtual void dg_ds(const State & state, SymSymR4 & res) const;
  /// Derivative of g wrt history
  virtual void dg_da(const State & state, History & res) const;

  /// Hardening rule proportional to the scalar strain rate
  virtual void h(const State & state, History & res) const;
  /// Derivative of h wrt stress
  virtual void dh_ds(const State & state, History & res) const;
  /// Derivative of h wrt history
  virtual void dh_da(const State & state, History & res) const;

 private:
  double eps0_, D_, n_, s0_, K_;
};

static Register<TestFlowRule> regTestFlowRule;

}
