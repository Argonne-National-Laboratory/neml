#include "singlecrystal.h"

namespace neml {

SingleCrystalModel::SingleCrystalModel(
    std::shared_ptr<KinematicModel> kinematics,
    std::shared_ptr<Lattice> lattice,
    std::shared_ptr<Orientation> initial_angle,
    std::shared_ptr<Interpolate> alpha,
    bool update_rotation, double tol, int miter, bool verbose, 
    int max_divide) :
      kinematics_(kinematics), lattice_(lattice), q0_(initial_angle), alpha_(alpha),
      update_rotation_(update_rotation), tol_(tol), miter_(miter),
      verbose_(verbose), max_divide_(max_divide), stored_hist_(false)
{
  populate_history(stored_hist_);
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
  pset.add_optional_parameter<int>("miter", 30);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<int>("max_divide", 6);

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
      params.get_parameter<bool>("verbose"),
      params.get_parameter<int>("max_divide"));
}

void SingleCrystalModel::populate_history(History & history) const
{
  history.add<Orientation>("rotation");
  history.add<Orientation>("rotation0");
  if (use_nye()) {
    history.add<RankTwo>("nye");
  }
  kinematics_->populate_history(history);
}

void SingleCrystalModel::init_history(History & history) const
{
  history.get<Orientation>("rotation") = *q0_;
  history.get<Orientation>("rotation0") = *q0_;
  if (use_nye()) {
    history.get<RankTwo>("nye") = RankTwo({0,0,0,0,0,0,0,0,0});
  }
  kinematics_->init_history(history);
}

double SingleCrystalModel::strength(double * const hist, double T) const
{
  History h_total = gather_history_(hist);
  History h = h_total.split(not_updated_());
  History f = h_total.split(not_updated_(), false);
  return kinematics_->strength(h, *lattice_, T, f);
}

void SingleCrystalModel::Fe(double * const stress, double * const hist,
                            double T, double * const Fe) const
{
  Symmetric S(stress);
  RankTwo FE(Fe);
  const History h = gather_history_(hist);

  Orientation Q = h.get<Orientation>("rotation").deepcopy();
  Orientation Q0 = h.get<Orientation>("rotation0").deepcopy();
 
  Orientation Re = Q * Q0.inverse();
  Symmetric estrain = kinematics_->elastic_strains(stress, Q, h, T);
  RankTwo RR;
  Re.to_matrix(RR.s());

  FE = ((Symmetric::id() + estrain) * RR).inverse();
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
  // Shut up a pointless memory error
  std::fill(h_np1, h_np1+nhist(), 0.0);

  // Setup everything in the appropriate wrappers
  const Symmetric D_np1(d_np1);
  const Skew W_np1(w_np1);
  const Symmetric D_n(d_n);
  const Skew W_n(w_n);

  double dt = t_np1 - t_n;
  
  // This is why I shouldn't do that this way
  Symmetric D;
  Skew W;
  if (dt != 0.0) {
    D = (D_np1 - D_n) / dt;
    W = (W_np1 - W_n) / dt;
  }

  Symmetric S_np1(s_np1);
  const Symmetric S_n(s_n);

  History HF_np1 = gather_history_(h_np1);
  const History HF_n = gather_history_(h_n);

  // There is a subtle race condition if you don't copy the lattice
  // because that object now caches data
  Lattice local_lattice = Lattice(*lattice_);

  // As the update is decoupled, split the histories into hardening/
  // orientation groups
  Orientation Q_n = HF_n.get<Orientation>("rotation");
  
  History H_np1 = HF_np1.split(not_updated_());
  History H_n = HF_n.split(not_updated_());
  History F_n = HF_n.split(not_updated_(), false);

  /* Begin adaptive stepping */
  double dT = T_np1 - T_n;
  int progress = 0;
  int cur_int_inc = pow(2, max_divide_);
  int target = cur_int_inc;
  int subdiv = 0;

  // Use S_np1 and H_np1 to iterate
  S_np1.copy_data(S_n.data());
  H_np1.copy_data(H_n.rawptr());
  
  while (progress < target) {
    double step = 1.0 / pow(2, subdiv);

    // Decouple the updates
    History fixed = kinematics_->decouple(S_np1, D, W, Q_n, H_np1, 
                                          local_lattice, T_n + dT *
                                          step, F_n);

    // Set the trial state
    SCTrialState trial(D, W,
                       S_np1, H_np1, // Yes, really
                       Q_n, local_lattice,
                       T_n + dT * step, dt * step,
                       fixed);

    // Solve the update
    int ier = solve_substep_(&trial, S_np1, H_np1);

    if (ier != 0) {
      subdiv++;
      cur_int_inc /= 2;

      if (verbose_) {
        std::cout << "Taking adaptive substep" << std::endl;
        std::cout << "Current progress " << progress << " out of " << target <<
            std::endl;
        std::cout << "New step fraction " << cur_int_inc << " out of "
            << target << std::endl;
      }

      if (subdiv >= max_divide_) {
        if (verbose_) {
          std::cout << "Adaptive substepping failed!" << std::endl;
        }
        return ier;
      }
    }
    else {
      progress += cur_int_inc;
      if (verbose_) {
        std::cout << "Adaptive substep succeeded" << std::endl;
        std::cout << "Current progress " << progress << " out of " << target <<
            std::endl;
      }
      
      // Calc tangent if we're going to be done
      if (progress == target) {
        // Tangent
        calc_tangents_(S_np1, H_np1, &trial, A_np1, B_np1);

        // Calculate the new rotation, if requested
        if (update_rotation_) {
          HF_np1.get<Orientation>("rotation") = update_rot_(S_np1, H_np1, &trial);
        }
        else {
          HF_np1.get<Orientation>("rotation") = Q_n;
        }
      }

    }
  }
  /* End adaptive stepping */

  /// If we store the Nye tensor just copy it over
  if (use_nye()) {
    //HF_np1.get<RankTwo>("nye") = HF_n.get<RankTwo>("nye");
  }

  // Calculate the new energy
  u_np1 = u_n + calc_energy_inc_(D_np1, D_n, S_np1, S_n);
  
  // Calculate the new dissipation
  p_np1 = p_n + calc_work_inc_(D_np1, D_n, S_np1, S_n, T_np1, T_n, 
                               HF_np1.get<Orientation>("rotation"), Q_n,
                               H_np1, H_n);

  return 0;
}

size_t SingleCrystalModel::nhist() const
{
  return stored_hist_.size();
}

int SingleCrystalModel::init_hist(double * const hist) const
{
  std::fill(hist, hist+nhist(), 0.0); // Shuts it up about initialized memory
  History h = gather_history_(hist);
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
  Symmetric stress(s_np1);

  const History h = gather_history_(h_np1);
  
  Symmetric estrain = kinematics_->elastic_strains(stress, 
                                                   h.get<Orientation>("rotation"), 
                                                   h, T_np1);
  std::copy(estrain.data(), estrain.data()+6, e_np1);

  return 0;
}

size_t SingleCrystalModel::nparams() const
{
  // 6 stress + the history minus the rotations minus the nye tensor
  if (use_nye()) {
    return 6 + nhist() - 8 - 9;
  }
  else {
    return 6 + nhist() - 8;
  }
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
  Symmetric S (x);
  History H = ats->history.copy_blank();
  H.copy_data(&x[6]);

  History & fixed = ats->fixed;

  // Get actual residual components
  Symmetric stress_res = S - ats->S - kinematics_->stress_rate(S, ats->d, ats->w, ats->Q,
                                                   H, ats->lattice, ats->T,
                                                   fixed) * ats->dt;
  History history_rate = kinematics_->history_rate(S, ats->d, ats->w, ats->Q,
                                                   H, ats->lattice, ats->T,
                                                   fixed);
  
  // Stick in residual
  std::copy(stress_res.data(), stress_res.data()+6, R);
  for (size_t i = 0; i < H.size(); i++) {
    R[i+6] = H.rawptr()[i] - ats->history.rawptr()[i] - history_rate.rawptr()[i] * ats->dt;
  }

  // Get all the Jacobian contributions
  SymSymR4 dSdS = kinematics_->d_stress_rate_d_stress(S, ats->d, ats->w, ats->Q,
                                                    H, ats->lattice, ats->T,
                                                    fixed);
  History dSdH = kinematics_->d_stress_rate_d_history(S, ats->d, ats->w, ats->Q,
                                                      H, ats->lattice, ats->T,
                                                      fixed);

  History dHdS = kinematics_->d_history_rate_d_stress(S, ats->d, ats->w, ats->Q,
                                                      H, ats->lattice, ats->T,
                                                      fixed);
  History dHdH = kinematics_->d_history_rate_d_history(S, ats->d, ats->w, ats->Q,
                                                       H, ats->lattice, ats->T,
                                                       fixed);

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
      J[CINDEX(i,(j+6),nparams())] = -dSdH.rawptr()[CINDEX(j,i,6)] * ats->dt;
    }
  }
  for (size_t i = 0; i<nh; i++) {
    for (size_t j = 0; j<6; j++) {
      J[CINDEX((i+6),j,nparams())] = -dHdS.rawptr()[CINDEX(i,j,6)] * ats->dt;
    }
  }
  for (size_t i = 0; i<nh; i++) {
    for (size_t j = 0; j<nh; j++) {
      J[CINDEX((i+6),(j+6),nparams())] = -dHdH.rawptr()[CINDEX(i,j,nh)] * ats->dt;
    }
  }

  // Add the identity
  for (size_t i = 0; i<nparams(); i++) {
    J[CINDEX(i,i,nparams())] += 1.0;
  }

  return 0;
}

Orientation SingleCrystalModel::get_active_orientation(
    double * const hist) const
{
  History h = gather_history_(hist);

  return get_active_orientation(h);
}

Orientation SingleCrystalModel::get_active_orientation(
    const History & hist) const
{
  return hist.get<Orientation>("rotation");
}

Orientation SingleCrystalModel::get_passive_orientation(
    double * const hist) const
{
  History h = gather_history_(hist);

  return get_passive_orientation(h);
}

Orientation SingleCrystalModel::get_passive_orientation(
    const History & hist) const
{
  return get_active_orientation(hist).inverse();
}

void SingleCrystalModel::set_active_orientation(
    double * const hist, const Orientation & q)
{
  History h = gather_history_(hist);

  set_active_orientation(h, q);
}

void SingleCrystalModel::set_active_orientation(
    History & hist, const Orientation & q)
{
  hist.get<Orientation>("rotation") = q;
  // This may need to be considered carefully, but for the cases of either
  // 1) Setting initial orientations through an external mechanism
  // 2) Recrystallization or twinning
  // resetting the reference angle is correct.
  hist.get<Orientation>("rotation0") = q;
}

void SingleCrystalModel::set_passive_orientation(
    double * const hist, const Orientation & q)
{
  History h = gather_history_(hist);

  set_passive_orientation(h, q);
}

void SingleCrystalModel::set_passive_orientation(
    History & hist, const Orientation & q)
{
  set_active_orientation(hist, q.inverse());
}

bool SingleCrystalModel::use_nye() const
{
  return kinematics_->use_nye();
}

void SingleCrystalModel::update_nye(double * const hist,
                                    const double * const nye) const
{
  if (use_nye()) {
    History h = gather_history_(hist);
    h.get<RankTwo>("nye") = RankTwo(std::vector<double>(nye,nye+9));
  }
}

History SingleCrystalModel::gather_history_(double * data) const
{
  History h = gather_blank_history_();
  h.set_data(data);
  return h;
}

History SingleCrystalModel::gather_history_(const double * data) const
{
  History h = gather_blank_history_();
  h.set_data(const_cast<double*>(data));
  return h;
}

History SingleCrystalModel::gather_blank_history_() const
{
  return stored_hist_;
}

void SingleCrystalModel::calc_tangents_(Symmetric & S, History & H,
                                        SCTrialState * ts,
                                        double * const A, double * const B)
{
  // Resetup x
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  std::copy(S.s(),S.s()+6,x);
  std::copy(H.rawptr(), H.rawptr()+H.size(),&x[6]);

  // Get the jacobian contributions
  double * R = new double[nparams()];
  double * J = new double[nparams()*nparams()];
  
  RJ(x, ts, R, J);

  size_t nh = nparams() - 6;

  if (nh > 0) {
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
        J12[CINDEX(i,j,nh)] = J[CINDEX(i,(j+6),nparams())];
      }
    }
   
    for (size_t i=0; i<nh; i++) {
      for (size_t j=0; j<6; j++) {
        J21[CINDEX(i,j,6)] = J[CINDEX((i+6),j,nparams())];
      }
    }

    for (size_t i=0; i<nh; i++) {
      for (size_t j=0; j<nh; j++) {
        J22[CINDEX(i,j,nh)] = J[CINDEX((i+6),(j+6),nparams())];
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
    History & fixed = ts->fixed;
    // Get the extra matrices
    SymSymR4 sd = kinematics_->d_stress_rate_d_d(S, ts->d, ts->w, ts->Q, 
                                               H, ts->lattice, ts->T, fixed) 
        + kinematics_->d_stress_rate_d_d_decouple(S, ts->d, ts->w,
                                                  ts->Q, H,
                                                  ts->lattice, ts->T, fixed);
    History hd = kinematics_->d_history_rate_d_d(S, ts->d, ts->w, ts->Q, 
                                               H, ts->lattice, ts->T, fixed);
    hd += kinematics_->d_history_rate_d_d_decouple(S, ts->d, ts->w,
                                                  ts->Q, H,
                                                  ts->lattice, ts->T, fixed);
    // hd is stored column major...
    double * hdt = new double[nh*6];
    for (size_t i = 0; i < nh; i++) {
      for (size_t j = 0; j < 6; j++) {
        hdt[CINDEX(i,j,6)] = hd.rawptr()[CINDEX(j,i,nh)];
      }
    }


    // Sadly need one more intermediate
    double * I1 = new double[36];
    mat_mat(6, 6, nh, M1, hdt, I1);
    delete [] hdt;
    for (size_t i=0; i<36; i++) {
      I1[i] = (sd.data()[i] - I1[i]); // Would mult by dt if things were sane
    }

    // Alright, do the final multiplication
    mat_mat(6,6,6, M2, I1, A);

    delete [] I1;

    // Do W
    SymSkewR4 sw = kinematics_->d_stress_rate_d_w(S, ts->d, ts->w, ts->Q, 
                                               H, ts->lattice, ts->T, fixed);
    sw += kinematics_->d_stress_rate_d_w_decouple(S, ts->d, ts->w,
                                                  ts->Q, H,
                                                  ts->lattice, ts->T, fixed);

    History hw = kinematics_->d_history_rate_d_w(S, ts->d, ts->w, ts->Q, 
                                               H, ts->lattice, ts->T, fixed);
    hw += kinematics_->d_history_rate_d_w_decouple(S, ts->d, ts->w,
                                                   ts->Q, H,
                                                   ts->lattice, ts->T, fixed);
    // This is transposed
    double * hwt = new double [nh*3];
    for (size_t i = 0; i < nh; i++) {
      for (size_t j = 0; j < 3; j++) {
        hwt[CINDEX(i,j,3)] = hw.rawptr()[CINDEX(j,i,nh)];
      }
    }

    // Sadly need one more intermediate
    double * I2 = new double[18];
    mat_mat(6, 3, nh, M1, hwt, I2);
    delete [] hwt;
    for (size_t i=0; i<18; i++) {
      I2[i] = (sw.data()[i] - I2[i]); // Would mult by dt if things were sane
    }

    // Alright, do the final multiplication
    mat_mat(6,3,6, M2, I2, B);

    delete [] I2;

    delete [] M1;
    delete [] M2;
  }
  // This is much simpler!
  else {
    delete [] R;

    invert_mat(J, 6);

    // Do D
    History & fixed = ts->fixed;
    // Get the extra matrices
    SymSymR4 sd = kinematics_->d_stress_rate_d_d(S, ts->d, ts->w, ts->Q, 
                                               H, ts->lattice, ts->T, fixed) 
        + kinematics_->d_stress_rate_d_d_decouple(S, ts->d, ts->w,
                                                  ts->Q, H,
                                                  ts->lattice, ts->T, fixed);
    
    // Multiply
    mat_mat(6,6,6, J, sd.data(), A);

    // Do W
    SymSkewR4 sw = kinematics_->d_stress_rate_d_w(S, ts->d, ts->w, ts->Q, 
                                               H, ts->lattice, ts->T, fixed);
    sw += kinematics_->d_stress_rate_d_w_decouple(S, ts->d, ts->w,
                                                  ts->Q, H,
                                                  ts->lattice, ts->T, fixed);

    // Multiply
    mat_mat(6,3,6,J,sw.data(), B);

    delete [] J;
  }
}

Orientation SingleCrystalModel::update_rot_(Symmetric & S_np1, History & H_np1,
                                            SCTrialState * ts) const
{
  Skew spin = kinematics_->spin(S_np1, ts->d, ts->w, ts->Q, H_np1, ts->lattice, ts->T,
                                ts->fixed);
  return wexp(spin * ts->dt) * ts->Q;
}

double SingleCrystalModel::calc_energy_inc_(
    const Symmetric & D_np1, const Symmetric & D_n, 
    const Symmetric & S_np1, const Symmetric & S_n) const
{
  return (S_np1 - S_n).contract(D_np1 - D_n) / 2.0;
}

double SingleCrystalModel::calc_work_inc_(
    const Symmetric & D_np1, const Symmetric & D_n,
    const Symmetric & S_np1, const Symmetric & S_n,
    double T_np1, double T_n, const Orientation & Q_np1,
    const Orientation & Q_n, const History & H_np1, 
    const History & H_n) const
{
  double dU = calc_energy_inc_(D_np1, D_n, S_np1, S_n);

  Symmetric e_np1 = kinematics_->elastic_strains(S_np1, Q_np1, H_np1, T_np1);
  Symmetric e_n = kinematics_->elastic_strains(S_n, Q_n, H_n, T_n);

  double dE = (S_np1 - S_n).contract(e_np1 - e_n) / 2.0;

  return dU - dE;
}

int SingleCrystalModel::solve_substep_(SCTrialState * ts,
                                       Symmetric & stress,
                                       History & hist)
{
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  int ier = solve(this, x, ts, tol_, miter_, verbose_);

  // Only dump into new stress and hist if we pass
  if (ier != 0) return ier;

  // Dump the results
  stress.copy_data(x);
  hist.copy_data(&x[6]);

  return ier;
}

std::vector<std::string> SingleCrystalModel::not_updated_() const
{
  if (use_nye()) {
    return {"rotation", "rotation0", "nye"};
  }
  else {
    return {"rotation", "rotation0"};
  }
}

} // namespace neml
