#include "singlecrystal.h"

namespace neml {

SingleCrystalModel::SingleCrystalModel(
    std::shared_ptr<KinematicModel> kinematics,
    std::shared_ptr<Lattice> lattice,
    std::shared_ptr<Orientation> initial_angle,
    std::shared_ptr<Interpolate> alpha,
    double tol, int miter, bool verbose) :
      kinematics_(kinematics), lattice_(lattice), q0_(initial_angle), alpha_(alpha),
      tol_(tol), miter_(miter), verbose_(verbose)
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
  pset.add_parameter<NEMLObject>("lattice");
  pset.add_optional_parameter<NEMLObject>("initial_rotation", 
                                          std::make_shared<Orientation>());
  pset.add_optional_parameter<NEMLObject>("alpha",
                                          std::make_shared<ConstantInterpolate>(0.0));
  pset.add_optional_parameter<double>("tol", 1.0e-6);
  pset.add_optional_parameter<int>("miter", 20);
  pset.add_optional_parameter<bool>("verbose", false);

  return pset;
}

std::unique_ptr<NEMLObject> SingleCrystalModel::initialize(ParameterSet & params)
{
  return neml::make_unique<SingleCrystalModel>(
      params.get_object_parameter<KinematicModel>("kinematics"),
      params.get_object_parameter<Lattice>("lattice"),
      params.get_object_parameter<Orientation>("initial_rotation"),
      params.get_object_parameter<Interpolate>("alpha"),
      params.get_parameter<double>("tol"),
      params.get_parameter<int>("miter"),
      params.get_parameter<bool>("verbose"));
}

void SingleCrystalModel::populate_history(History & history) const
{
  history.add<Orientation>("rotation");
  kinematics_->populate_history(history);
}

void SingleCrystalModel::init_history(History & history) const
{
  history.get<Orientation>("rotation") = *q0_;
  kinematics_->init_history(history);
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
  // Setup everything in the appropriate wrappers
  // Note that we need to fix the const thing
  Symmetric D_np1 = Symmetric(const_cast<double*>(d_np1));
  Skew W_np1      = Skew(const_cast<double*>(w_np1));

  Symmetric S_np1 = Symmetric(s_np1);
  Symmetric S_n = Symmetric(const_cast<double*>(s_n));

  History HF_np1 = gather_history_();
  HF_np1.set_data(h_np1);

  History HF_n = gather_history_();
  HF_n.set_data(const_cast<double*>(h_n));

  // As the update is decoupled, split the histories into hardening/
  // orientation groups
  Orientation Q_n = HF_n.get<Orientation>("rotation");
  History H_np1 = HF_np1.split({"rotation"});
  History H_n = H_n.split({"rotation"});
 
  // Decouple the updates
  kinematics_->decouple(S_np1, D_np1, W_np1, Q_n, H_n, *lattice_, T_np1);

  // Set the trial state
  SCTrialState trial(D_np1, W_np1, S_n, H_n, Q_n, *lattice_, T_np1, t_np1 - t_n);

  // Solve the update
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  int ier = solve(this, x, &trial, tol_, miter_, verbose_);
  if (ier != 0) return ier;

  // Dump the results
  std::copy(x, x+6, S_np1.s());
  H_np1.copy_data(&x[6]);

  // Calculate the tangents
  
  // Calculate the new rotation

  // Rotate through the tangents?
  
  // Calculate the new energy

  return 0;
}

size_t SingleCrystalModel::nhist() const
{
  History h = gather_history_();
  return h.size();
}

int SingleCrystalModel::init_hist(double * const hist) const
{
  History h = gather_history_();
  h.set_data(hist);
  init_history(h);
  return 0;
}

double SingleCrystalModel::alpha(double T) const
{
  return alpha_->value(T);
}

int SingleCrystalModel::elastic_strains(
    const double * const s_np1,
    double T_np1, const double * const h_np1,
    double * const e_np1) const
{
  Symmetric stress(const_cast<double*>(s_np1)); // Ugly need to fix

  History h = gather_history_();
  h.set_data(const_cast<double*>(h_np1));
  
  Symmetric estrain = kinematics_->elastic_strains(stress, 
                                                   h.get<Orientation>("rotation"), 
                                                   h, T_np1);
  std::copy(estrain.data(), estrain.data()+6, e_np1);

  return 0;
}

size_t SingleCrystalModel::nparams() const
{
  // 6 stress + the history minus the rotation
  return 6 + nhist() - 4;
}

int SingleCrystalModel::init_x(double * const x, TrialState * ts)
{
  SCTrialState * ats = static_cast<SCTrialState*>(ts);
  std::copy(ats->S.data(), ats->S.data()+6, x);
  std::copy(ats->history.rawptr(), ats->history.rawptr()+ats->history.size(), &x[6]);
  return 0;
}

int SingleCrystalModel::RJ(const double * const x, TrialState * ts,
                           double * const R, double * const J)
{
  // Cast trial state
  SCTrialState * ats = static_cast<SCTrialState*>(ts);

  // Make nice objects
  Symmetric S;
  std::copy(x,x+6,S.s());
  History H = ats->history.copy_blank();
  H.copy_data(&x[6]);

  // Get actual residual components
  Symmetric stress_res = S - ats->S - kinematics_->stress_rate(S, ats->d, ats->w, ats->Q,
                                                   H, ats->lattice, ats->T) * ats->dt;
  History history_rate = kinematics_->history_rate(S, ats->d, ats->w, ats->Q,
                                                   H, ats->lattice, ats->T);
  
  // Stick in residual
  std::copy(stress_res.data(), stress_res.data()+6, R);
  for (size_t i = 0; i < H.size(); i++) {
    R[i+6] = H.rawptr()[i] - ats->history.rawptr()[i] - history_rate.rawptr()[i] * ats->dt;
  }

  // Get all the Jacobian contributions
  SymSym dSdS = kinematics_->d_stress_rate_d_stress(S, ats->d, ats->w, ats->Q,
                                                    H, ats->lattice, ats->T);
  History dSdH = kinematics_->d_stress_rate_d_history(S, ats->d, ats->w, ats->Q,
                                                      H, ats->lattice, ats->T);

  History dHdS = kinematics_->d_history_rate_d_stress(S, ats->d, ats->w, ats->Q,
                                                      H, ats->lattice, ats->T);
  History dHdH = kinematics_->d_history_rate_d_history(S, ats->d, ats->w, ats->Q,
                                                       H, ats->lattice, ats->T);

  size_t nh = nparams() - 6;

  // Block them into the matrix
  for (size_t i = 0; i<6; i++) {
    for (size_t j = 0; j<6; j++) {
      J[CINDEX(i,j,nparams())] = -dSdS.data()[CINDEX(i,j,6)] * ats->dt;
    }
  }
  // This guy is stored column major
  for (size_t i = 0; i<6; i++) {
    for (size_t j = 0; j<nh; j++) {
      J[CINDEX(i,j+6,nparams())] = -dSdH.rawptr()[CINDEX(j,i,6)] * ats->dt;
    }
  }
  for (size_t i = 0; i<nh; i++) {
    for (size_t j = 0; j<6; j++) {
      J[CINDEX(i+6,j,nparams())] = -dHdS.rawptr()[CINDEX(i,j,6)] * ats->dt;
    }
  }
  // This guy is stored column major
  for (size_t i = 0; i<nh; i++) {
    for (size_t j = 0; j<nh; j++) {
      J[CINDEX(i+6,j+6,nparams())] = -dHdH.rawptr()[CINDEX(j,i,nh)] * ats->dt;
    }
  }

  // Add the identity
  for (size_t i = 0; i<nparams(); i++) {
    J[CINDEX(i,i,nparams())] += 1.0;
  }

  return 0;
}

History SingleCrystalModel::gather_history_() const
{
  History h(false);
  populate_history(h);
  return h;
}

} // namespace neml
