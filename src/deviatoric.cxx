#include "deviatoric.h"

namespace neml {

// Base class implementation
DeviatoricModel::DeviatoricModel()
{

}

DeviatoricModel::~DeviatoricModel()
{

}

// Linear elastic model
LEModel::LEModel()
{

}

LEModel::~LEModel()
{

}

size_t LEModel::nhist() const
{

  return 0;
}

int LEModel::update(
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{

  return 0;
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
    const double * const e_np1, const double * const e_n,
    double T_np1, double T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1) const
{

  return 0;
}


} // namespace neml
