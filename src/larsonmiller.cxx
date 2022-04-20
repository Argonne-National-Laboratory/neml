#include "larsonmiller.h"

#include "math/nemlmath.h"

namespace neml {

LarsonMillerRelation::LarsonMillerRelation(ParameterSet & params) :
    NEMLObject(params),
    fn_(params.get_object_parameter<Interpolate>("function")), 
    C_(params.get_parameter<double>("C")), 
    rtol_(params.get_parameter<double>("rtol")), 
    atol_(params.get_parameter<double>("atol")), 
    miter_(params.get_parameter<int>("miter")), 
    verbose_(params.get_parameter<bool>("verbose")),
    linesearch_(params.get_parameter<bool>("linesearch"))
{

}

std::string LarsonMillerRelation::type()
{
  return "LarsonMillerRelation";
}

ParameterSet LarsonMillerRelation::parameters()
{
  ParameterSet pset(LarsonMillerRelation::type());

  pset.add_parameter<NEMLObject>("function");
  pset.add_parameter<double>("C");
  pset.add_optional_parameter<double>("rtol", 1.0e-6);
  pset.add_optional_parameter<double>("atol", 1e-6);
  pset.add_optional_parameter<int>("miter", 20);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<bool>("linesearch", false);

  return pset;
}

std::unique_ptr<NEMLObject> LarsonMillerRelation::initialize(
    ParameterSet & params)
{
  return neml::make_unique<LarsonMillerRelation>(params); 
}

void LarsonMillerRelation::sR(double t, double T, double & s) const
{
  double LMP = T * (C_ + log10(t));
  s = pow(10.0,fn_->value(LMP));
}

void LarsonMillerRelation::tR(double s, double T, double & t)
{
  LMTrialState ts;
  ts.stress = s;
  double x[1];
  solve(this, x, &ts, {rtol_, atol_, miter_, verbose_, 
                  linesearch_});

  // x has LMP
  t = pow(10.0, x[0] / T - C_);
}

void LarsonMillerRelation::dtR_ds(double s, double T, double & dt)
{
  LMTrialState ts;
  ts.stress = s;
  double x[1];
  solve(this, x, &ts, {rtol_, atol_, miter_, verbose_, 
                  linesearch_});
  
  double LMP = x[0];
  
  dt = pow(10.0, LMP/T-C_)/(T * s * fn_->derivative(LMP));
}

size_t LarsonMillerRelation::nparams() const
{
  return 1;
}

void LarsonMillerRelation::init_x(double * const x, TrialState * ts)
{
  x[0] = 50000.0; // A reasonable LMP...
}

void LarsonMillerRelation::RJ(const double * const x, TrialState * ts, 
                             double * const R, double * const J)
{
  LMTrialState * tss = static_cast<LMTrialState*>(ts);

  R[0] = log10(tss->stress) - fn_->value(x[0]);
  J[0] = -fn_->derivative(x[0]);
}

}
