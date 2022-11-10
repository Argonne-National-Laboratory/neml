#pragma once

#include "../models.h"
#include "../math/rotations.h"
#include "singlecrystal.h"

#include "../windows.h"

namespace neml {

/// Generic superclass
class NEML_EXPORT PolycrystalModel: public NEMLModel_ldi
{
 public:
  PolycrystalModel(ParameterSet & params);

  size_t n() const;

  virtual void populate_hist(History & hist) const;
  virtual void init_hist(History & hist) const;

  virtual void populate_state(History & hist) const {};
  virtual void init_state(History & hist) const {};

  double * history(double * const store, size_t i) const;
  double * stress(double * const store, size_t i) const;
  double * d(double * const store, size_t i) const;
  double * w(double * const store, size_t i) const;

  const double * history(const double * const store, size_t i) const;
  const double * stress(const double * const store, size_t i) const;
  const double * d(const double * const store, size_t i) const;
  const double * w(const double * const store, size_t i) const;

  virtual std::vector<Orientation> orientations(double * const store) const;
  virtual std::vector<Orientation> orientations_active(double * const store) const;

 protected:
  std::shared_ptr<SingleCrystalModel> model_;
  const std::vector<std::shared_ptr<Orientation>> q0s_;
  int nthreads_;
  std::vector<double> weights_;
};

class NEML_EXPORT TaylorModel: public PolycrystalModel
{
 public:
  TaylorModel(ParameterSet & params);

  /// Type for the object system
  static std::string type();
  /// Parameters for the object system
  static ParameterSet parameters();
  /// Setup from a ParameterSet
  static std::unique_ptr<NEMLObject> initialize(ParameterSet & params);

  /// Large strain incremental update
  virtual void update_ld_inc(
     const double * const d_np1, const double * const d_n,
     const double * const w_np1, const double * const w_n,
     double T_np1, double T_n,
     double t_np1, double t_n,
     double * const s_np1, const double * const s_n,
     double * const h_np1, const double * const h_n,
     double * const A_np1, double * const B_np1,
     double & u_np1, double u_n,
     double & p_np1, double p_n);

  virtual double alpha(double T) const;
  virtual void elastic_strains(const double * const s_np1,
                              double T_np1, const double * const h_np1,
                              double * const e_np1) const;
};

static Register<TaylorModel> regTaylorModel;

}
