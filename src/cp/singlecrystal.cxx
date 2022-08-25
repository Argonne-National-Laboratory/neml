#include "cp/singlecrystal.h"

namespace neml {

SingleCrystalModel::SingleCrystalModel(ParameterSet & params) :
    NEMLModel_ldi(params),
    kinematics_(params.get_object_parameter<KinematicModel>("kinematics")), 
    lattice_(params.get_object_parameter<Lattice>("lattice")), 
    q0_(params.get_object_parameter<Orientation>("initial_rotation")), 
    alpha_(params.get_object_parameter<Interpolate>("alpha")),
    update_rotation_(params.get_parameter<bool>("update_rotation")), 
    rtol_(params.get_parameter<double>("rtol")), 
    atol_(params.get_parameter<double>("atol")), 
    miter_(params.get_parameter<int>("miter")),
    verbose_(params.get_parameter<bool>("verbose")),
    linesearch_(params.get_parameter<bool>("linesearch")),
    max_divide_(params.get_parameter<int>("max_divide")), 
    postprocessors_(params.get_object_parameter_vector<CrystalPostprocessor>("postprocessors")),
    elastic_predictor_(params.get_parameter<bool>("elastic_predictor")),
    fallback_elastic_predictor_(params.get_parameter<bool>("fallback_elastic_predictor")),
    force_divide_(params.get_parameter<int>("force_divide")),
    elastic_predictor_first_step_(params.get_parameter<bool>(
            "elastic_predictor_first_step"))
{
  cache_history_();
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
                                          zero_orientation());
  pset.add_optional_parameter<NEMLObject>("alpha",
                                          make_constant(0.0));
  pset.add_optional_parameter<bool>("update_rotation", true);
  pset.add_optional_parameter<double>("rtol", 1.0e-8);
  pset.add_optional_parameter<double>("atol", 1.0e-6);
  pset.add_optional_parameter<int>("miter", 30);
  pset.add_optional_parameter<bool>("verbose", false);
  pset.add_optional_parameter<bool>("linesearch", false);
  pset.add_optional_parameter<int>("max_divide", 6);
  pset.add_optional_parameter<std::vector<NEMLObject>>("postprocessors", {});
  pset.add_optional_parameter<bool>("elastic_predictor", false);
  pset.add_optional_parameter<bool>("fallback_elastic_predictor", true);
  pset.add_optional_parameter<int>("force_divide", 0);
  pset.add_optional_parameter<bool>("elastic_predictor_first_step", false);

  return pset;
}

std::unique_ptr<NEMLObject> SingleCrystalModel::initialize(ParameterSet & params)
{
  return neml::make_unique<SingleCrystalModel>(params);
}

void SingleCrystalModel::populate_state(History & history) const
{
  kinematics_->populate_hist(history);
}

void SingleCrystalModel::populate_static(History & history) const
{
  NEMLModel_ldi::populate_static(history);

  history.add<Orientation>("rotation");
  history.add<Orientation>("rotation0");
  if (use_nye()) {
    history.add<RankTwo>("nye");
  }
  for (auto pp : postprocessors_)
    pp->populate_hist(*lattice_, history);
}

void SingleCrystalModel::init_state(History & history) const
{
  kinematics_->init_hist(history);
}

void SingleCrystalModel::init_static(History & history) const
{
  NEMLModel_ldi::init_static(history);

  history.get<Orientation>("rotation") = *q0_;
  history.get<Orientation>("rotation0") = *q0_;
  if (use_nye()) {
    history.get<RankTwo>("nye") = RankTwo({0,0,0,0,0,0,0,0,0});
  }
  for (auto pp : postprocessors_)
    pp->init_hist(*lattice_, history);
}

double SingleCrystalModel::strength(double * const hist, double T) const
{
  History h_total = gather_history_(hist);
  auto [h, f] = split_state(h_total);
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
  Symmetric estrain = kinematics_->elastic_strains(stress, *lattice_, Q, h, T);
  RankTwo RR;
  Re.to_matrix(RR.s());

  FE = ((Symmetric::id() + estrain) * RR).inverse();
}

void SingleCrystalModel::update_ld_inc(
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
  // First step
  if ((t_n == 0.0) && elastic_predictor_first_step_) {
    attempt_update_ld_inc_(d_np1, d_n, w_np1, w_n,
                                    T_np1, T_n, t_np1, t_n,
                                    s_np1, s_n, h_np1, h_n, 
                                    A_np1, B_np1, u_np1, u_n,
                                    p_np1, p_n, 1);
    return;
  }

  // Try with a predictor
  if (elastic_predictor_) {
    attempt_update_ld_inc_(d_np1, d_n, w_np1, w_n,
                                    T_np1, T_n, t_np1, t_n,
                                    s_np1, s_n, h_np1, h_n, 
                                    A_np1, B_np1, u_np1, u_n,
                                    p_np1, p_n, 1);
  }
  else {
    // Base update (no elastic predictor)
    try {
      attempt_update_ld_inc_(d_np1, d_n, w_np1, w_n,
                                    T_np1, T_n, t_np1, t_n,
                                    s_np1, s_n, h_np1, h_n, 
                                    A_np1, B_np1, u_np1, u_n,
                                    p_np1, p_n, 0);
    }
    catch (const NEMLError & e) {
      // If it fails try with an elastic predictor
      if (fallback_elastic_predictor_)
        attempt_update_ld_inc_(d_np1, d_n, w_np1, w_n,
                               T_np1, T_n, t_np1, t_n,
                               s_np1, s_n, h_np1, h_n, 
                               A_np1, B_np1, u_np1, u_n,
                               p_np1, p_n, 1);
      else
        throw e;
    }
  }
}

void SingleCrystalModel::attempt_update_ld_inc_(
   const double * const d_np1, const double * const d_n,
   const double * const w_np1, const double * const w_n,
   double T_np1, double T_n,
   double t_np1, double t_n,
   double * const s_np1, const double * const s_n,
   double * const h_np1, const double * const h_n,
   double * const A_np1, double * const B_np1,
   double & u_np1, double u_n,
   double & p_np1, double p_n, int trial_type)
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
  
  auto [H_np1, F_np1] = split_state(HF_np1);
  auto [H_n, F_n] = split_state(HF_n);

  /* Begin adaptive stepping */
  double dT = T_np1 - T_n;
  int progress = 0;
  int cur_int_inc = pow(2, max_divide_-force_divide_);
  int target = pow(2, max_divide_);
  int subdiv = force_divide_;

  // Use S_np1 and H_np1 to iterate
  S_np1.copy_data(S_n.data());
  H_np1.copy_data(H_n.rawptr());
  
  while (progress < target) {
    double step = 1.0 / pow(2, subdiv);

    // Decouple the updates
    History fixed = kinematics_->decouple(S_np1, D, W, Q_n, H_np1, 
                                          local_lattice, T_n + dT *
                                          step, F_n);
    
    Symmetric strial;
    if (trial_type == 1) {
      // Predict the elastic unload
      strial = S_n + kinematics_->stress_increment(S_np1, D, W, dt * step,
                                                             local_lattice, Q_n,
                                                             H_np1,
                                                             T_n+dT*step);
    }
    else {
      strial = S_n;
    }
    
    // Set the trial state
    SCTrialState trial(D, W,
                       strial, S_n, H_np1, // Yes, really
                       Q_n, local_lattice,
                       T_n + dT * step, dt * step,
                       fixed);

    // Solve the update
    try {
      solve_substep_(&trial, S_np1, H_np1);
    }
    catch (const NEMLError & e) {
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
        throw NonlinearSolverError("Exceeded maximum adaptive subdivisions");
      }
      continue;
    }
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

  // Update model based on any post-processors
  for (auto pp : postprocessors_)
    pp->act(*this, local_lattice, T_np1, D, W, HF_np1, HF_n);
}

double SingleCrystalModel::alpha(double T) const
{
  return alpha_->value(T);
}

void SingleCrystalModel::elastic_strains(
    const double * const s_np1,
    double T_np1, const double * const h_np1,
    double * const e_np1) const
{
  Symmetric stress(s_np1);

  const History h = gather_history_(h_np1);
  
  Symmetric estrain = kinematics_->elastic_strains(stress, *lattice_,
                                                   h.get<Orientation>("rotation"), 
                                                   h, T_np1);
  std::copy(estrain.data(), estrain.data()+6, e_np1);
}

size_t SingleCrystalModel::nparams() const
{
  return 6 + nstate();
}

void SingleCrystalModel::init_x(double * const x, TrialState * ts)
{
  SCTrialState * ats = static_cast<SCTrialState*>(ts);
  std::copy(ats->S.data(), ats->S.data()+6, x);
  std::copy(ats->history.rawptr(), ats->history.rawptr()+ats->history.size(), &x[6]);
}

void SingleCrystalModel::RJ(const double * const x, TrialState * ts,
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
  Symmetric stress_res = S - ats->S_n - kinematics_->stress_rate(S, ats->d, ats->w, ats->Q,
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

  Symmetric e_np1 = kinematics_->elastic_strains(S_np1, *lattice_, Q_np1,
                                                 H_np1, T_np1);
  Symmetric e_n = kinematics_->elastic_strains(S_n, *lattice_, Q_n, H_n, T_n);

  double dE = (S_np1 - S_n).contract(e_np1 - e_n) / 2.0;

  return dU - dE;
}

void SingleCrystalModel::solve_substep_(SCTrialState * ts,
                                       Symmetric & stress,
                                       History & hist)
{
  std::vector<double> xv(nparams());
  double * x = &xv[0];
  solve(this, x, ts, {rtol_, atol_, miter_, verbose_, linesearch_});

  // Dump the results
  stress.copy_data(x);
  hist.copy_data(&x[6]);
}

} // namespace neml
