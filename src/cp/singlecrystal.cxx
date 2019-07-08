#include "singlecrystal.h"

namespace neml {

SingleCrystalModel::SingleCrystalModel(
    std::shared_ptr<KinematicModel> kinematics, 
    std::shared_ptr<Orientation> initial_angle,
    std::shared_ptr<Interpolate> alpha) :
      kinematics_(kinematics), q0_(initial_angle), alpha_(alpha)
{

}

SingleCrystalModel::~SingleCrystalModel()
{

}

std::string SingleCrystalModel::type()
{
  return "SingleCrystalModel";
}

ParameterSet SingleCrystalModel::parameters()
{
  ParameterSet pset(SingleCrystalModel::type());
  
  pset.add_parameter<NEMLObject>("kinematics");
  pset.add_optional_parameter<NEMLObject>("initial_rotation", 
                                          std::make_shared<Orientation>());
  pset.add_optional_parameter<NEMLObject>("alpha",
                                          std::make_shared<ConstantInterpolate>(0.0));

  return pset;
}

std::unique_ptr<NEMLObject> SingleCrystalModel::initialize(ParameterSet & params)
{
  return neml::make_unique<SingleCrystalModel>(
      params.get_object_parameter<KinematicModel>("kinematics"),
      params.get_object_parameter<Orientation>("initial_rotation"),
      params.get_object_parameter<Interpolate>("alpha"));
}

int SingleCrystalModel::update_ld_inc(
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
  return 0;
}

size_t SingleCrystalModel::nhist() const
{
  return 0;
}

int SingleCrystalModel::init_hist(double * const hist) const
{
  return 0;
}

double SingleCrystalModel::alpha(double T) const
{
  return 0.0;
}

int SingleCrystalModel::elastic_strains(
    const double * const s_np1,
    double T_np1, const double * const h_np1,
    double * const e_np1) const
{
  return 0;
}

size_t SingleCrystalModel::nparams() const
{
  return 0;
}

int SingleCrystalModel::init_x(double * const x, TrialState * ts)
{
  return 0;
}

int SingleCrystalModel::RJ(const double * const x, TrialState * ts,
                           double * const R, double * const J)
{
  return 0;
}

} // namespace neml
