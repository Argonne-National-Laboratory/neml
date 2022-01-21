#pragma once

#include "history.h"
#include "math/tensors.h"

namespace neml {

// A bunch of compile time junk to map derivative types
// I'm only going to cover combinations of Symmetric and double

template<class A, class B>
struct dtype;

template <>
struct dtype<double,double> {
  typedef double type;
};

template <>
struct dtype<double,Symmetric> {
  typedef Symmetric type;
};

template <>
struct dtype<Symmetric,double> {
  typedef Symmetric type;
};

template <>
struct dtype<Symmetric,Symmetric> {
  typedef SymSymR4 type;
};

/// Template class to handle internal variable evolution
template <class V>
class InternalVariable : public NEMLObject {
 public:
  // Black magic type mapping
  typedef typename dtype<V,V>::type VV;
  typedef typename dtype<V,double>::type Vs;
  typedef typename dtype<V,Symmetric>::type VS;

 public:
  // State object
  struct VariableState {
    V h;
    double a;
    double adot;
    double D;
    Symmetric s;
    Symmetric g;
    double T;
  };

 public:
  InternalVariable(ParameterSet & params);
  
  std::string name() const {return name_;};
  void set_name(std::string name) {name_ = name;};

  virtual V initial_value() = 0;

  virtual V ratep(VariableState & state) = 0;
  virtual VV d_ratep_d_h(VariableState & state) = 0;
  virtual Vs d_ratep_d_a(VariableState & state) = 0;
  virtual Vs d_ratep_d_adot(VariableState & state) = 0;
  virtual Vs d_ratep_d_D(VariableState & state) = 0;
  virtual VS d_ratep_d_s(VariableState & state) = 0;
  virtual VS d_ratep_d_g(VariableState & state) = 0;

  virtual V ratet(VariableState & state) = 0;
  virtual VV d_ratet_d_h(VariableState & state) = 0;
  virtual Vs d_ratet_d_a(VariableState & state) = 0;
  virtual Vs d_ratet_d_adot(VariableState & state) = 0;
  virtual Vs d_ratet_d_D(VariableState & state) = 0;
  virtual VS d_ratet_d_s(VariableState & state) = 0;
  virtual VS d_ratet_d_g(VariableState & state) = 0;

  virtual V rateT(VariableState & state) = 0;
  virtual VV d_rateT_d_h(VariableState & state) = 0;
  virtual Vs d_rateT_d_a(VariableState & state) = 0;
  virtual Vs d_rateT_d_adot(VariableState & state) = 0;
  virtual Vs d_rateT_d_D(VariableState & state) = 0;
  virtual VS d_rateT_d_s(VariableState & state) = 0;
  virtual VS d_rateT_d_g(VariableState & state) = 0;

 private:
  std::string name_;
};

template <class V>
InternalVariable<V>::InternalVariable(ParameterSet & params) :
    NEMLObject(params),
    name_(params.get_parameter<std::string>("name"))
{

}

/// Specialize template class definition for scalar internal variables
using ScalarInternalVariable = InternalVariable<double>;

/// Specialize template class definition for Symmetric internal variables
using SymmetricInternalVariable = InternalVariable<Symmetric>;

} // namespace neml
