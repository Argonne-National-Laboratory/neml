#include "polycrystal.h"

#include "batch.h"

namespace neml {

PolycrystalModel::PolycrystalModel(std::shared_ptr<SingleCrystalModel> model,
                                   std::vector<std::shared_ptr<Orientation>> qs,
                                   int nthreads) :
    model_(model), q0s_(qs), nthreads_(nthreads)
{

}

size_t PolycrystalModel::n() const
{
  return q0s_.size();
}

size_t PolycrystalModel::nhist() const
{
  return (model_->nstore() + 6 + 6 + 3) * n();
}

int PolycrystalModel::init_hist(double * const hist) const
{
  for (size_t i = 0; i < n(); i++) {
    model_->init_store(history(hist, i));
    model_->set_active_orientation(history(hist,i), *q0s_[i]);

    std::fill(stress(hist, i), stress(hist, i) + 6, 0);
    std::fill(d(hist, i), d(hist, i) + 6, 0);
    std::fill(w(hist, i), w(hist, i) + 3, 0);
  }

  return 0;
}

double * PolycrystalModel::history(double * const store, size_t i) const 
{
  return &(store[i*(model_->nstore())]);
}

double * PolycrystalModel::stress(double * const store, size_t i) const
{
  return &(store[n()*(model_->nstore()) + i * 6]);
}

double * PolycrystalModel::d(double * const store, size_t i) const
{
  return &(store[n()*(model_->nstore() + 6) + i * 6]);
}

double * PolycrystalModel::w(double * const store, size_t i) const
{
  return &(store[n()*(model_->nstore() + 6 + 6) + i * 3]); 
}

const double * PolycrystalModel::history(const double * const store, size_t i) const 
{
  return &(store[i*(model_->nstore())]);
}

const double * PolycrystalModel::stress(const double * const store, size_t i) const
{
  return &(store[n()*(model_->nstore()) + i * 6]);
}

const double * PolycrystalModel::d(const double * const store, size_t i) const
{
  return &(store[n()*(model_->nstore() + 6) + i * 6]);
}

const double * PolycrystalModel::w(const double * const store, size_t i) const
{
  return &(store[n()*(model_->nstore() + 6 + 6) + i * 3]); 
}

std::vector<Orientation> PolycrystalModel::orientations(double * const store) const
{
  std::vector<Orientation> res;
  
  get_orientation_passive_batch(*model_, n(), 
                                history(store, 0), res);
  return res;
}

TaylorModel::TaylorModel(std::shared_ptr<SingleCrystalModel> model,
                         std::vector<std::shared_ptr<Orientation>> qs,
                         int nthreads) :
    PolycrystalModel(model, qs, nthreads)
{

}

std::string TaylorModel::type()
{
  return "TaylorModel";
}

ParameterSet TaylorModel::parameters()
{
  ParameterSet pset(TaylorModel::type());
  
  pset.add_parameter<NEMLObject>("model");
  pset.add_parameter<std::vector<NEMLObject>>("qs");
  pset.add_optional_parameter<int>("nthreads", 1);

  return pset;
}

std::unique_ptr<NEMLObject> TaylorModel::initialize(ParameterSet & params)
{
  return neml::make_unique<TaylorModel>(
      params.get_object_parameter<SingleCrystalModel>("model"),
      params.get_object_parameter_vector<Orientation>("qs"),
      params.get_parameter<int>("nthreads"));
}

size_t TaylorModel::nstore() const
{
  return nhist();
}

int TaylorModel::init_store(double * const store) const
{
  return init_hist(store);
}

int TaylorModel::update_ld_inc(
   const double * const d_np1, const double * const d_n,
   const double * const w_np1, const double * const w_n,
   double T_np1, double T_n,
   double t_np1, double t_n,
   double * const s_np1, const double * const s_n,
   double * const h_np1, const double * const h_n,
   double * const A_np1, double * const B_np1,
   double & u_np1, double u_n,
   double & p_np1, double p_n)
{
  std::fill(s_np1, s_np1+6, 0);
  std::fill(A_np1, A_np1+36, 0);
  std::fill(B_np1, B_np1+18, 0);
  u_np1 = 0;
  p_np1 = 0;

  double * A_local = new double[36*n()];
  double * B_local = new double[18*n()];
  double * u_local = new double[n()];
  double * p_local = new double[n()];

  double * zero = new double[n()];
  std::fill(zero, zero+n(), 0.0);

  double * Ts_np1 = new double[n()];
  std::fill(Ts_np1, Ts_np1+n(), T_np1);
  double * Ts_n = new double[n()];
  std::fill(Ts_n, Ts_n+n(), T_n);

  for (size_t i = 0; i < n(); i++) {
    std::copy(d_np1, d_np1+6, d(h_np1, i));
    std::copy(w_np1, w_np1+3, w(h_np1, i));
  }

  int ier = evaluate_crystal_batch(*model_, n(), 
                         d(h_np1, 0), d(h_n, 0),
                         w(h_np1, 0), w(h_n, 0),
                         Ts_np1, Ts_n,
                         t_np1, t_n,
                         stress(h_np1, 0), stress(h_n, 0),
                         history(h_np1, 0), history(h_n, 0),
                         A_local, B_local,
                         u_local, zero,
                         p_local, zero, nthreads_);
  
  delete [] zero;
  delete [] Ts_np1;
  delete [] Ts_n;
                         
  for (size_t i = 0; i < n(); i++) {
    for (size_t j = 0; j < 6; j++) s_np1[j] += stress(h_np1, i)[j];
    for (size_t j = 0; j < 36; j++) A_np1[j] += A_local[i*36+j];
    for (size_t j = 0; j < 18; j++) B_np1[j] += B_local[i*18+j];
    u_np1 += u_local[i];
    p_np1 += p_local[i];
  }

  delete [] A_local;
  delete [] B_local;
  delete [] u_local;
  delete [] p_local;

  for (size_t j = 0; j < 6; j++) s_np1[j] /= n();
  for (size_t j = 0; j < 36; j++) A_np1[j] /= n();
  for (size_t j = 0; j < 18; j++) B_np1[j] /= n();

  u_np1 /= n();
  p_np1 /= n();

  u_np1 += u_n;
  p_np1 += p_n;

  return ier;
}

double TaylorModel::alpha(double T) const
{
  return model_->alpha(T);
}

int TaylorModel::elastic_strains(const double * const s_np1, 
                                 double T_np1, const double * const h_np1, 
                                 double * const e_np1) const
{
  std::fill(e_np1, e_np1+6, 0.0);
  double e_local[6];
  
  for (size_t i = 0; i < n(); i++) {
    model_->elastic_strains(stress(h_np1, i), T_np1, history(h_np1, i),
                            e_local);
    for (size_t j = 0; j < n(); j++) e_np1[j] += e_local[j];
  }

  for (size_t j = 0; j < n(); j++) e_np1[j] /= n();

  return 0;
}

}
