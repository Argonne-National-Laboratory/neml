#include "cp/batch.h"

#ifdef USE_OMP
#include <omp.h>
#endif

namespace neml {

void evaluate_crystal_batch(SingleCrystalModel & model, size_t n, 
                           const double * const d_np1, const double * const d_n, 
                           const double * const w_np1, const double * const w_n, 
                           const double * const T_np1, const double * const T_n, 
                           double t_np1, double t_n, 
                           double * const s_np1, const double * const s_n, 
                           double * const h_np1, const double * const h_n, 
                           double * const A_np1, double * const B_np1, 
                           double * const u_np1, const double * const u_n, 
                           double * const p_np1, const double * const p_n,
                           int nthreads)
{
  size_t nh = model.nstore();
  
  #ifdef USE_OMP
  omp_set_num_threads(nthreads);
  #endif

#ifdef USE_OMP
#pragma omp parallel for
#endif
  for (size_t i=0; i<n; i++) {
    model.update_ld_inc(&d_np1[i*6], &d_n[i*6], &w_np1[i*3], &w_n[i*3],
                                 T_np1[i], T_n[i], t_np1, t_n, &s_np1[i*6],
                                 &s_n[i*6], &h_np1[i*nh], &h_n[i*nh],
                                 &A_np1[i*36], &B_np1[i*18], 
                                 u_np1[i], u_n[i], p_np1[i], p_n[i]);
  }
}

void init_history_batch(SingleCrystalModel & model, size_t n, double * const hist)
{
  size_t nh = model.nstore();
  for (size_t i=0; i<n; i++) {
    model.init_store(&hist[i*nh]);
  }
}

void set_orientation_passive_batch(SingleCrystalModel & model, size_t n,
                                  double * const hist,
                                  std::vector<Orientation> orientations)
{
  if (orientations.size() != n) 
    throw NEMLError("Vector of orientations does not match batch size");
  
  size_t nh = model.nstore();

  for (size_t i=0; i<n; i++) {
    model.set_passive_orientation(&hist[i*nh], orientations[i]);
  }
}

void get_orientation_passive_batch(SingleCrystalModel & model, size_t n,
                                  double * const hist,
                                  std::vector<Orientation> & orientations)
{
  orientations.resize(n);
  size_t nh = model.nstore();

  for (size_t i=0; i<n; i++) {
    orientations[i] = model.get_passive_orientation(&hist[i*nh]);
  }
}

void get_orientation_active_batch(SingleCrystalModel & model, size_t n,
                                  double * const hist,
                                  std::vector<Orientation> & orientations)
{
  orientations.resize(n);
  size_t nh = model.nstore();

  for (size_t i=0; i<n; i++) {
    orientations[i] = model.get_active_orientation(&hist[i*nh]);
  }
}

} // namespace neml
