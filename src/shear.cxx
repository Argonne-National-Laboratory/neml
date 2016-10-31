#include "shear.h"

namespace neml {

ShearModulus::ShearModulus()
{

}

ShearModulus::~ShearModulus()
{

}


// Implementation of constant shear modulus model
ConstantShearModulus::ConstantShearModulus(double mu) :
    ShearModulus(), mu_(mu)
{

}

ConstantShearModulus::~ConstantShearModulus()
{

}

double ConstantShearModulus::modulus(double T_np1)
{
  return mu_;
}


} // namespace neml
