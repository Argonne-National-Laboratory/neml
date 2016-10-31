#include "deviatoric.h"

namespace neml {

// Base class implementation
DeviatoricModel::DeviatoricModel()
{

}

DeviatoricModel::~DeviatoricModel()
{

}




// Rate independent, associative flow model implementation
RIAFModel::RIAFModel()
{


}

RIAFModel::~RIAFModel()
{

}

size_t RIAFModel::nhist() const
{

  return 0;
}

int RIAFModel::update(
    const double* const e_inc,
    double T_np1, double T_inc,
    double t_np1, double t_inc,
    double * const h_np1, const double * const h_n,
    double * const s_np1, const double * const s_n,
    double * const A_np1) const
{

  return 0;
}


} // namespace neml
