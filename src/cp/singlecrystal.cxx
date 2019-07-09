#include "singlecrystal.h"

namespace neml {

SingleCrystalModel::SingleCrystalModel(
    std::shared_ptr<KinematicModel> kinematics,
    std::shared_ptr<Lattice> lattice,
    std::shared_ptr<Orientation> initial_angle,
    std::shared_ptr<Interpolate> alpha,
    bool update_rotation, double tol, int miter, bool verbose) :
      kinematics_(kinematics), lattice_(lattice), q0_(initial_angle), alpha_(alpha),
      update_rotation_(update_rotation), tol_(tol), miter_(miter),
      verbose_(verbose)
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
  pset.add_optional_parameter<bool>("update_rotation", true);
  pset.add_optional_parameter<double>("tol", 1.0e-6);
  pset.add_optional_parameter<int>("miter", 100);
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
      params.get_parameter<bool>("update_rotation"),
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
  Symmetric D_np1;
  D_np1.copy_data(d_np1);
  Skew W_np1;
  W_np1.copy_data(w_np1);

  Symmetric D_n;
  D_n.copy_data(d_n);
  Skew W_n;
  W_n.copy_data(w_n);

  double dt = t_np1 - t_n;

  // This is why I shouldn't do that this way
  Symmetric D;
  Skew W;
  if (dt != 0.0) {
    D = (D_np1 - D_n) / dt;
    W = (W_np1 - W_n) / dt;
  }

  Symmetric S_np1(s_np1);

  Symmetric S_n;
  S_n.copy_data(s_n);

  History HF_np1 = gather_history_();
  HF_np1.set_data(h_np1);

  History HF_n = HF_np1.copy_blank();
  HF_n.copy_data(h_n);

  // As the update is decoupled, split the histories into hardening/
  // orientation groups
  Orientation Q_n = HF_n.get<Orientation>("rotation");
  History H_np1 = HF_np1.split({"rotation"});
  History H_n = HF_n.split({"rotation"});

  // Decouple the updates
  kinematics_->decouple(S_n, D, W, Q_n, H_n, *lattice_, T_np1);

  // Set the trial state
  SCTrialState trial(D, W, S_n, H_n, Q_n, *lattice_, T_np1, dt);

  // Solve the update
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  int ier = solve(this, x, &trial, tol_, miter_, verbose_);
  if (ier != 0) return ier;

  // Dump the results
  S_np1.copy_data(x);
  H_np1.copy_data(&x[6]);

  // Calculate the tangents
  calc_tangents_(x, &trial, A_np1, B_np1);
  
  // Calculate the new rotation
  if (update_rotation_) {

  }

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
  Symmetric stress;
  stress.copy_data(s_np1);

  History h = gather_history_();
  h.copy_data(h_np1);
  
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
  S.copy_data(x);
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

void SingleCrystalModel::calc_tangents_(double * const x, SCTrialState * ts,
                                        double * const A, double * const B)
{
  // Get nice stress and history 
  Symmetric S;
  S.copy_data(x);
  History H = ts->history.copy_blank();
  H.copy_data(&x[6]);
  
  // Get the jacobian contributions
  double * R = new double[nparams()];
  double * J = new double[nparams()*nparams()];
  
  RJ(x, ts, R, J);

  size_t nh = nparams() - 6;

  double * J11 = new double[6*6];
  double * J12 = new double[6*nh];
  double * J21 = new double[nh*6];
  double * J22 = new double[nh*nh];

  for (size_t i=0; i<6; i++) {
    for (size_t j=0; j<6; j++) {
      J11[CINDEX(i,j,6)] = J[CINDEX(i,j,nparams())];
    }
  }

  for (size_t i=0; i<6; i++) {
    for (size_t j=0; j<nh; j++) {
      J12[CINDEX(i,j,nh)] = J[CINDEX(i,j+6,nparams())];
    }
  }
 
  for (size_t i=0; i<nh; i++) {
    for (size_t j=0; j<6; j++) {
      J21[CINDEX(i,j,6)] = J[CINDEX(i+6,j,nparams())];
    }
  }

  for (size_t i=0; i<nh; i++) {
    for (size_t j=0; j<nh; j++) {
      J22[CINDEX(i,j,nh)] = J[CINDEX(i+6,j+6,nparams())];
    }
  }

  delete [] R;
  delete [] J;

  // Let the games begin
  // J22 inverse
  invert_mat(J22, nh);

  // J12 J22_inv
  double * M1 = new double [6*nh];
  mat_mat(6, nh, nh, J12, J22, M1);

  // J12 J22_inv J21
  double * M2 = new double[6*6];
  mat_mat(6, 6, nh, M1, J21, M2);

  // J11 - J12 J22_inv J21
  for (size_t i = 0; i < 36; i++) {
    M2[i] = J11[i] - M2[i];
  }

  // Inverse of J11 - J12 J22_inv J21
  invert_mat(M2, 6);

  delete [] J11;
  delete [] J12;
  delete [] J21;
  delete [] J22;
  
  // Do D
  // Get the extra matrices
  SymSym sd = kinematics_->d_stress_rate_d_d(S, ts->d, ts->w, ts->Q, 
                                             H, ts->lattice, ts->T);
  History hd = kinematics_->d_history_rate_d_d(S, ts->d, ts->w, ts->Q, 
                                             H, ts->lattice, ts->T);

  // Sadly need one more intermediate
  double * I1 = new double[36];
  mat_mat(6, 6, nh, M1, hd.rawptr(), I1);
  for (size_t i=0; i<36; i++) {
    I1[i] = (sd.data()[i] - I1[i]); // Would mult by dt if things were sane
  }

  // Alright, do the final multiplication
  mat_mat(6,6,6, M2, I1, A);

  delete [] I1;

  // Do W
  SymSkew sw = kinematics_->d_stress_rate_d_w(S, ts->d, ts->w, ts->Q, 
                                             H, ts->lattice, ts->T);
  History hw = kinematics_->d_history_rate_d_w(S, ts->d, ts->w, ts->Q, 
                                             H, ts->lattice, ts->T);

  // Sadly need one more intermediate
  double * I2 = new double[18];
  mat_mat(6, 3, nh, M1, hw.rawptr(), I2);
  for (size_t i=0; i<18; i++) {
    I2[i] = (sw.data()[i] - I2[i]); // Would mult by dt if things were sane
  }

  // Alright, do the final multiplication
  mat_mat(6,3,6, M2, I2, B);

  delete [] I2;

  delete [] M1;
  delete [] M2;
}

} // namespace neml
