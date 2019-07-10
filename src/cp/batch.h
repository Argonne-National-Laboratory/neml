#ifndef BATCH_H
#define BATCH_H

#include "singlecrystal.h"

#include "../nemlerror.h"

namespace neml {

int evaluate_crystal_batch(SingleCrystalModel & model, size_t n, 
                           const double * const d_np1, const double * const d_n, 
                           const double * const w_np1, const double * const w_n, 
                           const double * const T_np1, const double * const T_n, 
                           double t_np1, double t_n, 
                           double * const s_np1, const double * const s_n, 
                           double * const h_np1, const double * const h_n, 
                           double * const A_np1, double * const B_np1, 
                           double * const u_np1, const double * const u_n, 
                           double * const p_np1, const double * const p_n,
                           int nthreads = 1);
int init_history_batch(SingleCrystalModel & model, size_t n, double * const hist);
int set_orientation_passive_batch(SingleCrystalModel & model, size_t n,
                                  double * const hist,
                                  std::vector<Orientation> orientations);
int get_orientation_passive_batch(SingleCrystalModel & model, size_t n,
                                  double * const hist,
                                  std::vector<Orientation> & orientations);

} // namespace neml

#endif // BATCH_H
