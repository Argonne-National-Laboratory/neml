#pragma once

#include "models.h"

namespace neml {

/// Block update in tensor notation
//  Update an entire block of models
//  Input data must be in row major order (i.e. nblock is the first axes)
//  Input and output must be as full tensors (not Mandel vectors)
NEML_EXPORT void block_evaluate(
    std::shared_ptr<NEMLModel> model,
    size_t nblock,
    const double * const e_np1, const double * const e_n,
    const double * const T_np1, const double * const T_n,
    double t_np1, double t_n,
    double * const s_np1, const double * const s_n,
    double * const h_np1, const double * const h_n,
    double * const A_np1,
    double * const u_np1, const double * const u_n,
    double * const p_np1, const double *  p_n);

/// Fast block tensor to Mandel converter
NEML_EXPORT void t2m(const double * const tensor, double * const mandel, size_t nblock);

/// Static data for t2m
NEML_EXPORT extern const double t2m_array[54];

/// Fast block Mandel to tensor converter
NEML_EXPORT void m2t(const double * const mandel, double * const tensor, size_t nblock);

/// Static data for m2t
NEML_EXPORT extern const double m2t_array[54];

/// Fast block rank 4 Mandel to rank 4 tensor convert
NEML_EXPORT void m42t4(const double * const mandel, double * const tensor, size_t nblock);

/// Static data for m42t4
NEML_EXPORT extern const double m42t4_array[2916];

} // namespace neml
