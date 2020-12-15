#include "larsonmiller.h"

#include "math/nemlmath.h"

namespace neml {

LarsonMillerRelation::LarsonMillerRelation(std::shared_ptr<Interpolate> fn, 
                                           double C, double rtol, double atol,
                                           int miter, bool verbose, 
                                           bool linesearch) :
    fn_(fn), C_(C), rtol_(rtol), atol_(atol), miter_(miter), verbose_(verbose),
    linesearch_(linesearch)
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
  return neml::make_unique<LarsonMillerRelation>(
      params.get_object_parameter<Interpolate>("function"),
      params.get_parameter<double>("C"),
      params.get_parameter<double>("rtol"),
      params.get_parameter<double>("atol"),
      params.get_parameter<int>("miter"),
      params.get_parameter<bool>("verbose"),
      params.get_parameter<bool>("linesearch")
      ); 
}

int LarsonMillerRelation::sR(double t, double T, double & s) const
{
  double LMP = T * (C_ + log10(t));
  s = pow(10.0,fn_->value(LMP));
  return 0;
}

int LarsonMillerRelation::tR(double s, double T, double & t)
{
  LMTrialState ts;
  ts.stress = s;
  double x[1];
  int ier = solve(this, x, &ts, {rtol_, atol_, miter_, verbose_, 
                  linesearch_});
  if (ier != 0) return ier;

  // x has LMP
  t = pow(10.0, x[0] / T - C_);

  return 0;
}

int LarsonMillerRelation::dtR_ds(double s, double T, double & dt)
{
  LMTrialState ts;
  ts.stress = s;
  double x[1];
  int ier = solve(this, x, &ts, {rtol_, atol_, miter_, verbose_, 
                  linesearch_});
  if (ier != 0) return ier;
  
  double LMP = x[0];
  
  dt = pow(10.0, LMP/T-C_)/(T * s * fn_->derivative(LMP));

  return 0;
}

size_t LarsonMillerRelation::nparams() const
{
  return 1;
}

int LarsonMillerRelation::init_x(double * const x, TrialState * ts)
{
  x[0] = 50000.0; // A reasonable LMP...

  return 0;
}

int LarsonMillerRelation::RJ(const double * const x, TrialState * ts, 
                             double * const R, double * const J)
{
  LMTrialState * tss = static_cast<LMTrialState*>(ts);

  R[0] = log10(tss->stress) - fn_->value(x[0]);
  J[0] = -fn_->derivative(x[0]);

  return 0;
}

}
