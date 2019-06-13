#include "kinematics.h"

namespace neml {

NoInelasticity::NoInelasticity()
{

}

NoInelasticity::~NoInelasticity()
{

}

std::string NoInelasticity::type()
{
  return "NoInelasticity";
}

ParameterSet NoInelasticity::parameters()
{
  ParameterSet pset(NoInelasticity::type());

  return pset;
}

std::unique_ptr<NEMLObject> NoInelasticity::initialize(ParameterSet & params)
{
  return neml::make_unique<NoInelasticity>(); 
}

Symmetric NoInelasticity::d_p(const Symmetric & stress, const Symmetric & d,
                              const Skew & w, const History & history,
                              const Lattice & lattice)
{
  return Symmetric();
}

Skew NoInelasticity::w_p(const Symmetric & stress, const Symmetric & d,
                         const Skew & w, const History & history,
                         const Lattice & lattice)
{
  return Skew();
}

} // namespace neml
