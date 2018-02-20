import numpy as np
import numpy.linalg as la

import scipy.interpolate as inter
import scipy.optimize as opt
from numpy.polynomial.legendre import leggauss

import numpy.random as ra

class Driver(object):
  """
    Superclass of all drivers, basically just sets up history and reports
    results.
  """
  def __init__(self, model, verbose = False, rtol = 1.0e-6, atol = 1.0e-10,
      miter = 25, T_init = 0.0, no_thermal_strain = False):
    """
      Parameters:
        model       material model to play with

      Optional:
        verbose     verbose output
        rtol        relative tolerance, where needed
        atol        absolute tolerance, where needed
        miter       maximum iterations, where needed
    """
    self.model = model
    self.verbose = verbose
    self.rtol = rtol
    self.atol = atol
    self.miter = miter
    self.nts = no_thermal_strain

    self.stress_int = [np.zeros((6,))]

    self.stored_int = [self.model.init_store()]
    self.T_int = [T_init]
    self.t_int = [0.0]

    self.u_int = [0.0]
    self.p_int = [0.0]

  @property
  def stress(self):
    return np.array(self.stress_int)

  @property
  def stored(self):
    return np.array(self.stored_int)

  @property
  def history(self):
    return self.stored[:,:self.model.nhist]

  @property
  def T(self):
    return np.array(self.T_int)

  @property
  def t(self):
    return np.array(self.t_int)

  @property
  def u(self):
    return np.array(self.u_int)

  @property
  def p(self):
    return np.array(self.p_int)


class Driver_sd(Driver):
  """
    Superclass of generic small strain drivers, contains generic step methods.
  """
  def __init__(self, *args, **kwargs):
    super(Driver_sd, self).__init__(*args, **kwargs)
    self.strain_int = [np.zeros((6,))]
    self.thermal_strain_int = [np.zeros((6,))]
    self.mechanical_strain_int = [np.zeros((6,))]

  def solve_try(self, RJ, x0, extra = []):
    def s1(x0i):
      try:
        x = newton(RJ, x0i, verbose = self.verbose,
            rtol = self.rtol, atol = self.atol, miter = self.miter)
        return x, True
      except Exception:
        return np.zeros((12,)), False

    def s2(x0i):
      try:
        res = opt.root(RJ, x0i, jac = True, method = 'lm')
        return res.x, res.success
      except Exception:
        return np.zeros((12,)), False

    def s3(x0i):
      try:
        x = newton(RJ, x0i, verbose = self.verbose,
            rtol = self.rtol, atol = self.atol, miter = self.miter,
            linesearch = 'backtracking')
        return x, True
      except Exception:
        return np.zeros((12,)), False

    solvers = [s1,s2,s3]
    guesses = [x0] + extra
    
    success = False
    for xi in guesses:
      for solv in solvers:
        x, success = solv(xi)
        if success:
          break
      if success:
        break

    if not success:
      raise MaximumIterations()

    return x

  @property
  def strain(self):
    return np.array(self.strain_int)

  @property
  def thermal_strain(self):
    return np.array(self.thermal_strain_int)

  @property
  def mechanical_strain(self):
    return np.array(self.mechanical_strain_int)

  def update_thermal_strain(self, T_np1):
    if self.nts:
      return np.zeros((6,))
    else:
      dT = T_np1 - self.T_int[-1]
      a_np1 = self.model.alpha(T_np1)
      a_n = self.model.alpha(self.T_int[-1])

      return self.thermal_strain_int[-1] + dT * (a_np1 + a_n) / 2 * np.array([1.0,1,1,0,0,0])

  def strain_step(self, e_np1, t_np1, T_np1):
    """
      Take a strain-controlled step

      Parameters:
        e_np1       next strain
        t_np1       next time
        T_np1       next temperature
    """
    enext = self.update_thermal_strain(T_np1)
    s_np1, h_np1, A_np1, u_np1, p_np1 = self.model.update_sd(e_np1 - enext, 
        self.mechanical_strain_int[-1],
        T_np1, self.T_int[-1], t_np1, self.t_int[-1], self.stress_int[-1],
        self.stored_int[-1], self.u_int[-1], self.p_int[-1])

    self.strain_int.append(np.copy(e_np1))
    self.mechanical_strain_int.append(e_np1 - enext)
    self.thermal_strain_int.append(enext)
    self.stress_int.append(np.copy(s_np1))
    self.stored_int.append(np.copy(h_np1))
    self.T_int.append(T_np1)
    self.t_int.append(t_np1)
    self.u_int.append(u_np1)
    self.p_int.append(p_np1)

  def stress_step(self, s_np1, t_np1, T_np1):
    """
      Take a stress-controlled step

      Parameters:
        s_np1       next stress
        t_np1       next time
        T_np1       next temperature
    """
    enext = self.update_thermal_strain(T_np1)
    def RJ(e):
      s, h, A, u, p = self.model.update_sd(e - enext, self.mechanical_strain_int[-1],
        T_np1, self.T_int[-1], t_np1, self.t_int[-1], 
        self.stress_int[-1],
        self.stored_int[-1], self.u_int[-1], self.p_int[-1])
      R = s - s_np1
      return R, A

    #e_np1 = newton(RJ, self.strain_int[-1], verbose = self.verbose,
    #    rtol = self.rtol, atol = self.atol, miter = self.miter)
    if len(self.strain_int) > 1:
      inc = self.strain_int[-1] - self.strain_int[-2]
      extra = [self.strain_int[-1] + inc]
    else:
      extra = []

    e_np1 = self.solve_try(RJ, self.strain_int[-1], extra = extra)

    self.strain_step(e_np1, t_np1, T_np1)

  def erate_step(self, sdir, erate, t_np1, T_np1,
      einc_guess = None, ainc_guess = None):
    """
      Drive in a given stress direction at a prescribed strain rate, like
      an actual "stress controlled" experiment.

      Parameters:
        sdir        stress direction
        erate       strain rate (in the direction)
        t_np1       next time
        T_np1       next temperature

      Optional:
        einc_guess      a guess at the strain increment
        ainc_guess      a guess at the stress increment
    """
    sdir = sdir / la.norm(sdir)
    dt = t_np1 - self.t_int[-1]
    enext = self.update_thermal_strain(T_np1)

    def RJ(x):
      a = x[0]
      e_inc = x[1:]
      s, h, A, u, p = self.model.update_sd(self.strain_int[-1] + e_inc - enext,
          self.mechanical_strain_int[-1],
          T_np1, self.T_int[-1], t_np1, self.t_int[-1], self.stress_int[-1],
          self.stored_int[-1],
          self.u_int[-1], self.p_int[-1])

      R = np.zeros((7,))
      J = np.zeros((7,7))
      R[:6] = s - (sdir * a + self.stress_int[-1])
      R[6] = np.dot(e_inc, sdir) / dt - erate

      J[:6,0] = -sdir
      J[:6,1:] = A
      J[6,0] = 0.0
      J[6,1:] = sdir / dt

      return R, J
    
    x0 = np.zeros((7,))
    
    if einc_guess is not None:
      x0[1:] = einc_guess
    else:
      x0[1:] = sdir / 10000.0

    if ainc_guess is not None:
      x0[0] = ainc_guess
    else:
      x0[0] = 1.0

    #x = newton(RJ, x0, verbose = self.verbose,
    #    rtol = self.rtol, atol = self.atol, miter = self.miter)
    x = self.solve_try(RJ, x0)
    e_np1 = self.strain_int[-1] + x[1:]

    self.strain_step(e_np1, t_np1, T_np1)

    return x[1:], x[0]

  def erate_einc_step(self, sdir, erate, einc, T_np1, **kwargs):
    """
      Similar to erate_step but specify the strain increment instead of the
      time increment.

      Parameters:
        sdir        stress direction
        erate       strain rate, in stress direction
        einc        strain increment, in stress direction
        T_np1       temperature at next time step
    """
    dt = einc / erate
    return self.erate_step(sdir, erate, self.t_int[-1] + dt, T_np1, **kwargs)
  
  def srate_sinc_step(self, sdir, srate, sinc, T_np1):
    """
      Similar to rate_step but specify the stress increment instead of the
      time increment.

      Parameters:
        sdir        stress direction
        srate       stress rate
        sinc        stress increment
        T_np1       temperature at next time step

    """
    if np.allclose(sdir, 0.0):
      s_np1 = self.stress_int[-1]
    else:
      s_np1 = self.stress_int[-1] + sdir / la.norm(sdir) * sinc
    if np.isclose(srate, 0.0):
      dt = 0.0
    else:
      dt = np.abs(np.dot(s_np1 - self.stress_int[-1], sdir) / srate)
    
    self.stress_step(s_np1, self.t_int[-1] + dt, T_np1)

  def strain_hold_step(self, i, t_np1, T_np1):
    """
      A special, mixed step which holds the strain in index i constant
      while holding the stress in the other directions to their previous
      values

      Parameters:
        i           index to hold
        t_np1       next time
        T_np1       next temperature
    """
    enext = self.update_thermal_strain(T_np1)
    oset = sorted(list(set(range(6)) - set([i])))
    def RJ(e_np1):
      s, h, A, u, p = self.model.update_sd(e_np1 - enext,
          self.mechanical_strain_int[-1],
          T_np1, self.T_int[-1], t_np1, self.t_int[-1], self.stress_int[-1],
          self.stored_int[-1], self.u_int[-1], self.p_int[-1])

      R = np.zeros((6,))
      R[0] = e_np1[i] - self.strain_int[-1][i]
      R[1:] = s[oset] - self.stress_int[-1][oset]

      J = np.zeros((6,6))
      J[0,0] = 1.0
      J[1:,1:] = A[oset,:][:,oset]

      return R, J
    
    x0 = np.copy(self.strain_int[-1])
    #e_np1 = newton(RJ, x0, verbose = self.verbose,
    #    rtol = self.rtol, atol = self.atol, miter = self.miter)
    e_np1 = self.solve_try(RJ, x0)

    self.strain_step(e_np1, t_np1, T_np1)

class Driver_sd_twobar(Driver_sd):
  def __init__(self, A1, A2, T1, T2, period, P, load_time, *args, **kwargs):
    """
      Generic twobar driver.

      Parameters:
        A1          area of bar 1
        A2          area of bar 2
        T1          function giving temperature versus time for 1 cycle (bar 1)
        T2          function giving temperature versus time for 1 cycle (bar 2)
        period      actual period of a cycle
        P           load applied
        load_time   how much time to use to load

        model       material model to use

      Optional:
        nsteps_load     number of steps to load up in
        nsteps_cycle    number of steps per cycle
    """
    nsteps_load = kwargs.pop('nsteps_load', 50)
    nsteps_cycle = kwargs.pop('nsteps_cycle', 201)
    dstrain = kwargs.pop('dstrain', lambda t: 0.0)
    dts = kwargs.pop('dts', None)

    if (T1(0) != T2(0)):
      raise ValueError("Initial bar temperatures do not match!")

    super(Driver_sd_twobar, self).__init__(*args, 
        no_thermal_strain = True, **kwargs)

    self.A1 = A1
    self.A2 = A2
    self.T1_fn = T1
    self.T2_fn = T2
    self.period = period
    self.P = P

    self.load_time = load_time

    self.nsteps_load = nsteps_load
    self.nsteps_cycle = nsteps_cycle
    self.dstrain = dstrain

    self.strain1_int = [np.zeros((6,))]
    self.strain2_int = [np.zeros((6,))]
    self.strain1_mech_int = [np.zeros((6,))]
    self.strain2_mech_int = [np.zeros((6,))]
    
    self.strain1_plastic_int = [np.zeros((6,))]
    self.strain2_plastic_int = [np.zeros((6,))]
    
    self.strain1_thermal_int = [np.zeros((6,))]
    self.strain2_thermal_int = [np.zeros((6,))]

    self.stress1_int = [np.zeros((6,))]
    self.stress2_int = [np.zeros((6,))]

    self.store1_int = [self.model.init_store()]
    self.store2_int = [self.model.init_store()]

    self.u1_int = [0.0]
    self.u2_int = [0.0]

    self.p1_int = [0.0]
    self.p2_int = [0.0]

    self.T1_int = [T1(0.0)]
    self.T2_int = [T2(0.0)]

    self.T1_0 = T1(0.0)
    self.T2_0 = T2(0.0)

    self.dts = dts

  @property
  def T1(self):
    return np.array(self.T1_int)

  @property
  def T2(self):
    return np.array(self.T2_int)

  @property
  def strain1(self):
    return np.array(self.strain1_int)

  @property
  def strain2(self):
    return np.array(self.strain2_int)

  @property
  def strain1_mech(self):
    return np.array(self.strain1_mech_int)

  @property
  def strain2_mech(self):
    return np.array(self.strain2_mech_int)

  @property
  def strain1_plastic(self):
    return np.array(self.strain1_plastic_int)

  @property
  def strain2_plastic(self):
    return np.array(self.strain2_plastic_int)

  @property
  def stress1(self):
    return np.array(self.stress1_int)

  @property
  def stress2(self):
    return np.array(self.stress2_int)

  def take_step(self, P, T1, T2, dt, dstrain):
    """
      Take a single load step to get to P, T1, and T2 from the previous state.
    """
    def RJ(x):
      e1 = x[:6]
      e2 = x[6:]

      a1 = self.model.alpha(T1)
      a2 = self.model.alpha(T2)

      a1p = self.model.alpha(self.T1_int[-1])
      a2p = self.model.alpha(self.T2_int[-1])

      dT1 = T1 - self.T1_int[-1]
      dT2 = T2 - self.T2_int[-1]

      eT1 = self.strain1_thermal_int[-1] + dT1 * (a1+a1p)/2 * np.array([1,1,1,0,0,0]) 
      eT2 = self.strain2_thermal_int[-1] + dT2 * (a2+a2p)/2 * np.array([1,1,1,0,0,0])
      
      em1 = e1 - eT1
      em2 = e2 - eT2
      
      s1_np1, h1_np1, A1_np1, u1_np1, p1_np1 = self.model.update_sd(
          em1, self.strain1_mech_int[-1],
          T1, self.T1_int[-1], dt + self.t_int[-1], self.t_int[-1],
          self.stress1_int[-1], self.store1_int[-1], self.u1_int[-1], 
          self.p1_int[-1])
      s2_np1, h2_np1, A2_np1, u2_np1, p2_np1 = self.model.update_sd(
          em2, self.strain2_mech_int[-1],
          T2, self.T2_int[-1], dt + self.t_int[-1], self.t_int[-1],
          self.stress2_int[-1], self.store2_int[-1], self.u2_int[-1], 
          self.p2_int[-1])

      R = np.zeros((12,))
      R[0] = s1_np1[0] * self.A1 + s2_np1[0] * self.A2 - P
      R[1] = e1[0] - e2[0] - dstrain
      R[2:7] = s1_np1[1:]
      R[7:] = s2_np1[1:]

      J = np.zeros((12,12))

      J[0,:6] = A1_np1[0,:] * self.A1 
      J[0,6:] = A2_np1[0,:] * self.A2

      J[1,:6] = np.array([1,0,0,0,0,0])
      J[1,6:] = np.array([-1,0,0,0,0,0])

      J[2:7,:6] = A1_np1[1:,:]
      J[2:7,6:] = 0.0

      J[7:,:6] = 0.0
      J[7:,6:] = A2_np1[1:,:]
      
      return R, J
    
    x0 = np.hstack((self.strain1_int[-1], self.strain2_int[-1]))
    if len(self.strain1_int) > 1:
      xn = np.hstack((self.strain1_int[-2], self.strain2_int[-2])) 
    else:
      xn = np.zeros((12,))

    xguesses = []
    xguesses.append(np.copy(x0))
    xguesses.append(x0 + (x0 - xn))
    xguesses.append(np.zeros((12,)))
    xguesses.append(x0 - (x0 - xn))
    
    solvers = []
    def s1(x0i):
      try:
        x = newton(RJ, x0i, verbose = self.verbose,
            rtol = self.rtol, atol = self.atol, miter = self.miter)
        return x, True
      except Exception:
        return np.zeros((12,)), False

    def s5(x0i):
      try:
        x = newton(RJ, x0i, verbose = self.verbose,
            rtol = self.rtol, atol = self.atol, miter = self.miter,
            linesearch = 'backtracking')
        return x, True
      except Exception:
        return np.zeros((12,)), False

    def s2(x0i):
      try:
        res = opt.root(RJ, x0i, jac = True, method = 'lm')
        return res.x, res.success
      except Exception:
        return np.zeros((12,)), False

    def s3(x0i):
      try:
        res = opt.root(RJ, x0i, jac = True, method = 'hybr')
        return res.x, res.success
      except Exception:
        return np.zeros((12,)), False

    def s4(x0i):
      try:
        x = quasi_newton(RJ, x0i, verbose = self.verbose,
            rtol = self.rtol, atol = self.atol, miter = self.miter)
        return x, True
      except Exception:
        return np.zeros((12,)), False

    solvers = [s1, s5, s2]

    for i,xi in enumerate(xguesses):
      bme = False
      for j,sv in enumerate(solvers):
        x, success = sv(xi)
        if success:
          bme = True
          break
      if bme:
        break

    if not success:
      raise MaximumIterations()
    
    e1 = x[:6]
    e2 = x[6:]

    self.strain1_int.append(e1)
    self.strain2_int.append(e2)

    a1 = self.model.alpha(T1)
    a2 = self.model.alpha(T2)

    a1p = self.model.alpha(self.T1_int[-1])
    a2p = self.model.alpha(self.T2_int[-1])

    dT1 = T1 - self.T1_int[-1]
    dT2 = T2 - self.T2_int[-1]

    eT1 = self.strain1_thermal_int[-1] + dT1 * (a1+a1p)/2 * np.array([1,1,1,0,0,0]) 
    eT2 = self.strain2_thermal_int[-1] + dT2 * (a2+a2p)/2 * np.array([1,1,1,0,0,0])

    em1 = e1 - eT1
    em2 = e2 - eT2

    s1_np1, h1_np1, A1_np1, u1_np1, p1_np1 = self.model.update_sd(
        em1, self.strain1_mech_int[-1],
        T1, self.T1_int[-1], dt + self.t_int[-1], self.t_int[-1],
        self.stress1_int[-1], self.store1_int[-1], self.u1_int[-1], 
        self.p1_int[-1])
    s2_np1, h2_np1, A2_np1, u2_np1, p2_np1 = self.model.update_sd(
        em2, self.strain2_mech_int[-1],
        T2, self.T2_int[-1], dt + self.t_int[-1], self.t_int[-1],
        self.stress2_int[-1], self.store2_int[-1], self.u2_int[-1], 
        self.p2_int[-1])

    self.strain1_mech_int.append(em1)
    self.strain2_mech_int.append(em2)

    self.strain1_thermal_int.append(eT1)
    self.strain2_thermal_int.append(eT2)

    self.t_int.append(self.t_int[-1] + dt)
    self.T1_int.append(T1)
    self.T2_int.append(T2)

    self.u1_int.append(u1_np1)
    self.u2_int.append(u2_np1)
    self.u_int.append(u1_np1*self.A1 + u2_np1*self.A2)

    self.p1_int.append(p1_np1)
    self.p2_int.append(p2_np1)
    self.p_int.append(p1_np1*self.A1 + p2_np1*self.A2)

    self.store1_int.append(h1_np1)
    self.store2_int.append(h2_np1)

    self.stress1_int.append(s1_np1)
    self.stress2_int.append(s2_np1)

    estrain1 = self.model.elastic_strains(s1_np1, T1, h1_np1)
    estrain2 = self.model.elastic_strains(s2_np1, T2, h2_np1)

    self.strain1_plastic_int.append(em1 - estrain1)
    self.strain2_plastic_int.append(em2 - estrain2)

  def load_up(self):
    """
      Load up to stress
    """
    dt = self.load_time / self.nsteps_load
    for P in np.linspace(0, self.P, self.nsteps_load+1)[1:]:
      self.take_step(P, self.T1_0, self.T2_0, dt, 0)

  def one_cycle(self):
    """
      Do one cycle at fixed load
    """
    if self.dts is None:
      dt = self.period / self.nsteps_cycle
      dts = [dt] * self.nsteps_cycle
    else:
      dts = self.dts
    
    t = 0.0
    for dt in dts:
      t += dt
      self.take_step(self.P, self.T1_fn(t), self.T2_fn(t), dt,
          self.dstrain(t))

def twobar_test(model, A1, A2, T1, T2, period, P, load_time, ncycles,
    max_strain = None, min_strain = None, 
    nsteps_load = 50, nsteps_cycle = 201, verbose = False,
    rtol_classify = 1.0e-4, atol_classify = 1.0e-10,
    dstrain = lambda t: 0.0, dts = None, ret_all = False):
  """
    Run a two bar test and classify

    Parameters:
      model         material model
      A1            area of bar 1
      A2            area of bar 2
      T1            one cycle temperature history, bar 1
      T2            one cycle temperature history, bar 2
      period        one cycle's time
      P             applied load
      load_time     physical time of loading
      ncycles       number of cycles to run

    Optional:
      max_strain    an absolute inelastic strain value to stop cycling at
      nsteps_load   number of steps to take to get to load
      nsteps_cycle  number of steps to take to take in a cycle
      verbose       print extra messages
      rtol_classify relative tolerance for cycle classification
      atol_classify absolute tolerance for cycle classification
      dstrain       direct strain difference between bars
      dts           manually specified cycle time steps
      ret_all       don't exclude the loading parts of the curves
  """
  if dts is not None:
    nsteps_cycle = len(dts)

  driver = Driver_sd_twobar(A1, A2, T1, T2, period, P, load_time, model,
      nsteps_load = nsteps_load, nsteps_cycle = nsteps_cycle, verbose = verbose,
      dstrain = dstrain, dts = dts)
  
  try:
    driver.load_up()
  except Exception:
    return {
        'strain': [],
        'ncycles': 0,
        'T1': [],
        'T2': [],
        'time': [],
        'energy_density': [], 
        'plastic_work': [],
        'strain_mech_1': [],
        'strain_mech_2': [],
        'strain_inelastic_1': [],
        'strain_inelastic_2': [],
        'stress1': [],
        'stress2': [],
        'classification': "collapse",
        'strain1': [],
        'strain2': []
        }
  
  cycle_count = 0
  unstable = False
  for i in range(ncycles):
    try:
      driver.one_cycle()
    except Exception as err:
      unstable = True
      break

    ep1 = driver.strain1_plastic_int[-1][0]
    ep2 = driver.strain2_plastic_int[-1][0]
    max_s = np.max([np.abs(ep1), np.abs(ep2)])
    min_s = np.min([np.abs(ep1), np.abs(ep2)])
    cycle_count += 1
    if (max_strain is not None) and (max_s > max_strain):
      break
    if (min_strain is not None) and (min_s > min_strain):
      break
  
  if unstable:
    classification = 'unstable'
  else:
    # Use the final cycles
    a = cycle_count * nsteps_cycle - 1
    b = (cycle_count - 1) * nsteps_cycle - 1

    classification = classify(
        driver.u_int[a], driver.u_int[b],
        driver.p_int[a], driver.p_int[b],
        driver.strain1_int[a][0], driver.strain1_int[b][0],
        driver.strain2_int[a][0], driver.strain2_int[b][0],
        rtol = rtol_classify,
        atol = atol_classify)

  if ret_all:
    nsteps_load = 0
  
  return {
      'strain': driver.strain1[nsteps_load:], 
      'ncycles': i+1,
      'T1': driver.T1[nsteps_load:],
      'T2': driver.T2[nsteps_load:],
      'time': driver.t[nsteps_load:] - driver.t[nsteps_load],
      'energy_density': driver.u[nsteps_load:], 
      'plastic_work': driver.p[nsteps_load:],
      'strain_mech_1': driver.strain1_mech[nsteps_load:,0],
      'strain_mech_2': driver.strain2_mech[nsteps_load:,0],
      'strain_inelastic_1': driver.strain1_plastic[nsteps_load:,0],
      'strain_inelastic_2': driver.strain2_plastic[nsteps_load:,0],
      'stress1': driver.stress1[nsteps_load:,0],
      'stress2': driver.stress2[nsteps_load:,0],
      'classification': classification,
      'strain1': driver.strain1[nsteps_load:,0],
      'strain2': driver.strain2[nsteps_load:,0],
      'nsteps_cycle': nsteps_cycle,
      'energy': driver.u_int[nsteps_load:],
      'dissipation': driver.p_int[nsteps_load:]
      }

def classify(ua, ub, pa, pb, e1a, e1b, e2a, e2b, rtol = 1.0e-4, atol = 1.0e-10):
  """
    Classify a model as elastic, elastic shakedown, plastic shakedown,
    or ratcheting.

    Parameters:
      ua    cycle a internal energy
      ub    cycle b internal energy
      pa    cycle a plastic dissipation
      pb    cycle b plastic dissipation
      ea    cycle a strain
      eb    cycle b strain

    Optional:
      rtol  relative tolerance
      atol  absolute tolerance
  """
  if np.abs(pb) < atol:
    return 'elastic'
  elif np.abs(ub-ua) < rtol * np.abs(ub):
    return 'elastic shakedown'
  elif (np.abs(e1b-e1a) < rtol * np.abs(e1b)) and (np.abs(e2b-e2a) < rtol * np.abs(e2b)):
    return 'plastic shakedown'
  else:
    return 'ratcheting'

def uniaxial_test(model, erate, T = 300.0, emax = 0.05, nsteps = 250, 
    sdir = np.array([1,0,0,0,0,0]), verbose = False,
    offset = 0.2/100.0, history = None):
  """
    Make a uniaxial stress/strain curve

    Parameters:
      model     material model
      erate     strain rate
    
    Optional:
      T         temperature, default 300.0
      emax      maximum strain, default 5%
      nsteps    number of steps to use, default 250
      sdir      stress direction, default tension in x
      verbose   whether to be verbose
      offset    used to calculate yield stress
      history   initial model history

    Results:
      strain    strain in direction
      stress    stress in direction
  """
  e_inc = emax / nsteps
  driver = Driver_sd(model, verbose = verbose, T_init = T)
  if history is not None:
    driver.stored_int[0] = history

  strain = [0.0]
  stress = [0.0]
  for i in range(nsteps):
    if i == 0:
      einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T)
    else:
      einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T, 
          einc_guess = einc, ainc_guess = ainc)
    strain.append(np.dot(driver.strain_int[-1], sdir))
    stress.append(np.dot(driver.stress_int[-1], sdir))

  strain = np.array(strain)
  stress = np.array(stress)

  # Calculate the yield stress and Young's modulus
  E = np.abs(stress[1]) / np.abs(strain[1])
  sfn = inter.interp1d(np.abs(strain), np.abs(stress))
  tfn = lambda e: E * (e - offset)
  try:
    sYe = opt.brentq(lambda e: sfn(e) - tfn(e), 0.0, np.max(strain))
    sY = tfn(sYe)
  except Exception:
    sY = np.inf

  return {'strain': strain, 'stress': stress, 
      'energy_density': np.copy(driver.u),
      'plastic_work': np.copy(driver.p),
      'youngs': E, 'yield': sY}

def strain_cyclic(model, emax, R, erate, ncycles, T = 300.0, nsteps = 50,
    sdir = np.array([1,0,0,0,0,0]), hold_time = None, n_hold = 25,
    verbose = False, check_dmg = False, dtol = 0.75):
  """
    Strain controlled cyclic test.

    Parameters:
      emax      maximum strain
      R         R = emin / emax
      erate     strain rate to go at
      ncycles   number of cycles

    Optional:
      T         temperature, default 300
      nsteps    number of steps per half cycle
      sdir      stress direction, defaults to x and tension first
      hold_time if None don't hold, if scalar then hold symmetrically top/bot
                if an array specify different hold times for first direction
                (default tension) and second direction
      n_hold    number of steps to hold over
      verbose   whether to be verbose

    Results:
      strain        strain in direction
      stress        stress in direction
      cycles        list of cycle numbers
      max           maximum stress per cycle
      min           minimum stress per cycle
      mean          mean stress per cycle
      
  """
  # Setup
  driver = Driver_sd(model, verbose = verbose, T_init = T)
  emin = emax * R
  if hold_time:
    if np.isscalar(hold_time):
      hold_time = [hold_time, hold_time]
  else:
    hold_time = [0,0]

  # Setup results
  strain = [0.0]
  stress = [0.0]
  time = [0.0]
  cycles = []
  smax = []
  smin = []
  smean = []

  ecycle = []
  pcycle = []

  # First half cycle
  if verbose:
    print("Initial half cycle")
  e_inc = emax / nsteps
  try:
    for i in range(nsteps):
      if i == 0:
        einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T)
      else:
        einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T, einc_guess = einc,
            ainc_guess = ainc)
      if check_dmg:
        if driver.stored_int[-1][0] > dtol:
          raise Exception("Damage check exceeded")
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))
      time.append(time[-1] + e_inc / erate)
  except Exception as e:
    print("Failed to make first half cycle")
    raise e
  
  # Begin cycling
  for s in range(ncycles):
    if verbose:
      print("Cycle %i" % s)

    try:
      # Tension hold
      if hold_time[0] > 0.0:
        dt = hold_time[0] / n_hold
        for i in range(n_hold):
          einc, ainc = driver.erate_step(sdir, 0.0, time[-1] + dt, T, 
              einc_guess = np.zeros((6,)), ainc_guess = -1)
          if check_dmg:
            if driver.stored_int[-1][0] > dtol:
              raise Exception("Damage check exceeded")
          strain.append(np.dot(driver.strain_int[-1], sdir))
          stress.append(np.dot(driver.stress_int[-1], sdir))
          time.append(time[-1] + dt)

      si = len(driver.strain_int)
      e_inc = np.abs(emin - emax) / nsteps
      for i in range(nsteps):
        if i == 0:
          einc, ainc = driver.erate_einc_step(-sdir, erate, e_inc, T, 
              einc_guess = -einc, ainc_guess = -ainc)
        else:
          einc, ainc = driver.erate_einc_step(-sdir, erate, e_inc, T, 
              einc_guess = einc, ainc_guess = ainc)
        if check_dmg:
          if driver.stored_int[-1][0] > dtol:
            raise Exception("Damage check exceeded")
        strain.append(np.dot(driver.strain_int[-1], sdir))
        stress.append(np.dot(driver.stress_int[-1], sdir))
        time.append(time[-1] + e_inc / erate)
      
      # Compression hold
      if hold_time[1] > 0.0:
        dt = hold_time[1] / n_hold
        for i in range(n_hold):
          einc, ainc = driver.erate_step(sdir, 0.0, time[-1] + dt, T, 
              einc_guess = np.zeros((6,)), ainc_guess = -1)
          if check_dmg:
            if driver.stored_int[-1][0] > dtol:
              raise Exception("Damage check exceeded")
          strain.append(np.dot(driver.strain_int[-1], sdir))
          stress.append(np.dot(driver.stress_int[-1], sdir))
          time.append(time[-1] + dt)

      e_inc = np.abs(emax - emin) / nsteps
      for i in range(nsteps):
        if i == 0:
          einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T, 
              einc_guess = -einc, ainc_guess = -ainc)
        else:
          einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T,
              einc_guess = einc, ainc_guess = ainc)
        if check_dmg:
          if driver.stored_int[-1][0] > dtol:
            raise Exception("Damage check exceeded")
        strain.append(np.dot(driver.strain_int[-1], sdir))
        stress.append(np.dot(driver.stress_int[-1], sdir))
        time.append(time[-1] + e_inc / erate)

      # Calculate
      if np.isnan(max(stress[si:])) or np.isnan(min(stress[si:])):
        break
      
      cycles.append(s)
      smax.append(max(stress[si:]))
      smin.append(min(stress[si:]))
      smean.append((smax[-1]+smin[-1])/2)

      ecycle.append(driver.u_int[-1])
      pcycle.append(driver.p_int[-1])
    except Exception as e:
      break

  # Setup and return
  return {"strain": np.array(strain), "stress": np.array(stress),
      "cycles": np.array(cycles, dtype = int), "max": np.array(smax),
      "min": np.array(smin), "mean": np.array(smean),
      "energy_density": np.array(ecycle), "plastic_work": np.array(pcycle),
      "history": driver.stored_int[-1], "time": np.array(time)}

def stress_cyclic(model, smax, R, srate, ncycles, T = 300.0, nsteps = 50,
    sdir = np.array([1,0,0,0,0,0]), hold_time = None, n_hold = 10,
    verbose = False, etol = 0.1):
  """
    Stress controlled cyclic test.

    Parameters:
      smax      maximum stress
      R         R = smin / smax
      srate     stress rate to go at
      ncycles   number of cycles

    Optional:
      T         temperature, default 300
      nsteps    number of steps per half cycle
      sdir      stress direction, defaults to x and tension first
      hold_time if None don't hold, if scalar then hold symmetrically top/bot
                if an array specify different hold times for first direction
                (default tension) and second direction
      n_hold    number of steps to hold over
      verbose   whether to be verbose

    Results:
      strain        strain in direction
      stress        stress in direction
      cycles        list of cycle numbers
      max           maximum strain per cycle
      min           minimum strain per cycle
      mean          mean stress per cycle
      
  """
  # Setup
  driver = Driver_sd(model, verbose = verbose, T_init = T)
  smin = smax * R
  if hold_time:
    if np.isscalar(hold_time):
      hold_time = [hold_time, hold_time]

  # Setup results
  strain = [0.0]
  stress = [0.0]
  cycles = []
  emax = []
  emin = []
  emean = []

  ecycle = []
  pcycle = []

  # First half cycle
  s_inc = smax / nsteps
  if verbose:
    print("First half cycle")
  for i in range(nsteps):
    driver.srate_sinc_step(sdir, srate, s_inc, T)
    strain.append(np.dot(driver.strain_int[-1], sdir))
    stress.append(np.dot(driver.stress_int[-1], sdir))

  # Begin cycling
  for s in range(ncycles):
    quit = False
    if verbose:
      print("Cycle %i" % s)
    si = len(driver.strain_int)
    # Hold, if requested
    if hold_time and (hold_time[0] > 0.0) and (s != 0):
      ht = hold_time[0]
      dt = ht / n_hold
      for i in range(n_hold):
        try:
          driver.stress_step(driver.stress_int[-1], driver.t_int[-1] + dt, T)
        except:
          quit = True
          break
        if la.norm(driver.strain_int[-1] - driver.strain_int[-2]) > etol:
          quit = True
          break
        strain.append(np.dot(driver.strain_int[-1], sdir))
        stress.append(np.dot(driver.stress_int[-1], sdir))
    
    if quit:
      break
    
    s_inc = (smin - smax) / nsteps
    for i in range(nsteps):
      try:
        driver.srate_sinc_step(sdir, srate, s_inc, T)
      except:
        quit = True
        break
      if la.norm(driver.strain_int[-1] - driver.strain_int[-2]) > etol:
        quit = True
        break
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

    if quit:
      break

    # Hold, if requested
    if hold_time and (hold_time[1] > 0.0):
      ht = hold_time[1]
      dt = ht / n_hold
      for i in range(n_hold):
        try:
          driver.stress_step(driver.stress_int[-1], driver.t_int[-1] + dt, T)
        except:
          quit = True
          break
        if la.norm(driver.strain_int[-1] - driver.strain_int[-2]) > etol:
          quit = True
          break
        strain.append(np.dot(driver.strain_int[-1], sdir))
        stress.append(np.dot(driver.stress_int[-1], sdir))

    if quit:
      break

    s_inc = (smax - smin) / nsteps
    for i in range(nsteps):
      try:
        driver.srate_sinc_step(sdir, srate, s_inc, T)
      except:
        quit = True
        break
      if la.norm(driver.strain_int[-1] - driver.strain_int[-2]) > etol:
        quit = True
        break
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

    if quit:
      break

    # Calculate
    cycles.append(s)
    emax.append(max(strain[si:]))
    emin.append(min(strain[si:]))
    emean.append((emax[-1]+emin[-1])/2)
    ecycle.append(driver.u_int[-1])
    pcycle.append(driver.p_int[-1])

  # Setup and return
  return {"strain": np.array(strain), "stress": np.array(stress),
      "cycles": np.array(cycles, dtype = int), "max": np.array(emax),
      "min": np.array(emin), "mean": np.array(emean),
      "energy_density": np.array(ecycle), "plastic_work": np.array(pcycle),
      "time": np.array(driver.t_int)}

def stress_relaxation(model, emax, erate, hold, T = 300.0, nsteps = 750,
    nsteps_up = 150, index = 0, tc = 1.0,
    verbose = False):
  """
    Simulate a stress relaxation test.

    Parameters:
      model         material model
      emax          maximum strain to attain
      erate         strain rate to take getting there
      hold          hold time

    Optional:
      T             temperature
      nsteps        number of steps to relax over
      nsteps_up     number of steps to take getting up to stress
      index         direction to pull in, default x tension
      tc            1.0 for tension -1.0 for compression
      verbose       whether to be verbose

    Results:
      time          time
      strain        strain
      stress        stress
      rtime         relaxation time
      rrate         stress relaxation rate
  """
  # Setup
  driver = Driver_sd(model, verbose = verbose, T_init = T)
  time = [0]
  strain = [0]
  stress = [0]

  # Ramp up
  if verbose:
    print("Ramp up")
  sdir = np.zeros((6,))
  sdir[index] = tc
  einc = emax / nsteps_up
  for i in range(nsteps_up):
    if i == 0:
      eincg, ainc = driver.erate_einc_step(sdir, erate, einc, T)
    else:
      eincg, ainc = driver.erate_einc_step(sdir, erate, einc, T,
          einc_guess = eincg, ainc_guess = ainc)
    time.append(driver.t[-1])
    strain.append(np.dot(driver.strain_int[-1],sdir))
    stress.append(np.dot(driver.stress_int[-1],sdir))

  ri = len(driver.strain_int)
  
  if verbose:
    print("Hold")
  dt = hold / nsteps
  for i in range(nsteps):
    driver.strain_hold_step(index, driver.t_int[-1] + dt, T)
    time.append(driver.t_int[-1])
    strain.append(np.dot(driver.strain_int[-1],sdir))
    stress.append(np.dot(driver.stress_int[-1],sdir))

  time = np.array(time)
  strain = np.array(strain)
  stress = np.array(stress)
  rrate = -np.diff(stress[ri:]) / np.diff(time[ri:])
  
  return {'time': np.copy(time), 'strain': np.copy(strain), 
      'stress': np.copy(stress), 'rtime': np.copy(time[ri:-1] - time[ri]),
      'rrate': np.copy(rrate), 'rstress': np.copy(stress[ri:-1])}


def creep(model, smax, srate, hold, T = 300.0, nsteps = 250,
    nsteps_up = 150, sdir = np.array([1,0,0,0,0,0]), verbose = False,
    logspace = False, history = None, elimit = 1.0, check_dmg = False,
    dtol = 0.75):
  """
    Simulate a creep test

    Parameters:
      model         material model
      smax          stress to attain
      srate         stress rate to take getting there
      hold          total hold time

    Optional:
      T             temperature
      nsteps        number of steps over relaxation period
      nsteps_up     number of steps to get to stress value
      sdir          stress direction, defaults to x-tension
      verbose       whether to be verbose
      logspace      if true logspace the time steps
      history       use damaged material
      check_dmg     check damage as a break condition
      dtol          damage to define failure at
  """
  # Setup
  driver = Driver_sd(model, verbose = verbose, T_init = T)
  if history is not None:
    driver.stored_int[0] = history
  time = [0]
  strain = [0]
  stress = [0]

  # Ramp up
  sinc = float(smax) / nsteps_up
  for i in range(nsteps_up):
    driver.srate_sinc_step(sdir, srate, sinc, T)
    time.append(driver.t[-1])
    strain.append(np.dot(driver.strain_int[-1],sdir))
    stress.append(np.dot(driver.stress_int[-1],sdir))

  ri = len(driver.strain_int)
  
  t0 = time[-1]
  if logspace:
    ts = np.logspace(0, np.log10(hold), num = nsteps) + t0
  else:
    ts = np.linspace(0,hold, num = nsteps) + t0
    
  failed = False
  for t in ts:
    # You can exceed the creep life of the sample doing this...
    # Need to allow a non-convergent result to break
    try:
      driver.stress_step(driver.stress_int[-1], t, T)
    except:
      failed = True
      break
    if np.any(np.isnan(driver.strain_int[-1])):
      failed = True
      break
    if np.any(np.abs(driver.strain_int[-1]) > elimit):
      failed = True
      break
    
    ed = np.dot(driver.strain_int[-1],sdir)
    if ed < strain[-1]:
      failed = True
      break

    if check_dmg:
      if driver.stored_int[-1][0] > dtol:
        failed = True
        break

    time.append(t)
    strain.append(ed)
    stress.append(np.dot(driver.stress_int[-1],sdir))

  time = np.array(time)
  strain = np.array(strain)
  stress = np.array(stress)
  rrate = np.diff(strain[ri:]) / np.diff(time[ri:])
  if len(strain) > ri +1:
    rstrain = strain[ri:] - strain[ri]
    rtime = time[ri:] - time[ri]
  else:
    rstrain = []
    rtime = []
  
  return {'time': np.copy(time), 'strain': np.copy(strain), 
      'stress': np.copy(stress), 'rtime': np.copy(rtime[:-1]),
      'rrate': np.copy(rrate), 'rstrain': np.copy(rstrain[:-1]),
      'tstrain': np.copy(strain[ri:-1]), 'failed': failed}

def thermomechanical_strain_raw(model, time, temperature, strain, 
    sdir = np.array([1,0,0,0,0,0.0]), verbose = False, substep = 1):
  """
    Directly drive a model using the output of a strain controlled thermomechanical test

    Parameters:
      model         material model
      time          list of times
      temperature   list of temperatures
      strain        list of strains

    Optional:
      sdir          direction of stress
      verbose       print information
      substep       take substep steps per data point
  """
  stress = np.zeros((len(time),))
  mechstrain = np.zeros((len(time),))
  driver = Driver_sd(model, verbose = verbose, T_init = temperature[0])

  einc = None
  ainc = None

  for i in range(1,len(stress)):
    quit = False
    for k in range(substep):
      ei_np1 = (strain[i] - strain[i-1]) / substep * (k+1) + strain[i-1]
      ti_np1 = (time[i] - time[i-1]) / substep * (k+1) + time[i-1]
      Ti_np1 = (temperature[i] - temperature[i-1]) / substep * (k+1) + temperature[i-1]

      ei_n = (strain[i] - strain[i-1]) / substep * (k) + strain[i-1]
      ti_n = (time[i] - time[i-1]) / substep * (k) + time[i-1]
      Ti_n = (temperature[i] - temperature[i-1]) / substep * (k) + temperature[i-1]

      erate = (ei_np1 - ei_n) / (ti_np1 - ti_n)
      try:
        if i == 1:
          einc, ainc = driver.erate_step(sdir, erate, ti_np1, Ti_np1)
        else:
          einc, ainc = driver.erate_step(sdir, erate, ti_np1, Ti_np1,
              einc_guess = einc, ainc_guess = ainc)
      except MaximumIterations:
        quit = True
        break
    if quit:
      break
    
    stress[i] = np.dot(driver.stress_int[-1], sdir)
    mechstrain[i] = np.dot(driver.thermal_strain_int[-1], sdir)


  return {'time': np.copy(time)[:i], 'temperature': np.copy(temperature)[:i], 'strain': np.copy(strain)[:i],
      'stress': np.copy(stress)[:i], 'mechanical strain': np.copy(mechstrain)[:i]}


def rate_jump_test(model, erates, T = 300.0, e_per = 0.01, nsteps_per = 100, 
    sdir = np.array([1,0,0,0,0,0]), verbose = False, history = None, 
    strains = None):
  """
    Model a uniaxial strain rate jump test

    Parameters:
      model         material model
      erate         list of strain rates
    
    Optional:
      T             temperature, default 300.0
      e_per         how much straining to do for each rate
      nsteps_per    number of steps per strain rate
      sdir          stress direction, default tension in x
      verbose       whether to be verbose
      history       prior model history
      strains       manual definition of jumps

    Results:
      strain    strain in direction
      stress    stress in direction
  """
  e_inc = e_per / nsteps_per
  driver = Driver_sd(model, verbose = verbose, T_init = T)
  if history is not None:
    driver.stored_int[0] = history
  strain = [0.0]
  stress = [0.0]
  
  if strains is None:
    for erate in erates:
      for i in range(nsteps_per):
        if i == 0:
          einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T)
        else:
          einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T, 
              einc_guess = einc, ainc_guess = ainc)
        strain.append(np.dot(driver.strain_int[-1], sdir))
        stress.append(np.dot(driver.stress_int[-1], sdir))
  else:
    incs = np.diff(np.insert(strains,0,0)) / nsteps_per
    for e,erate,inc in zip(strains,erates,incs):
      while strain[-1] < e:
        einc, ainc = driver.erate_einc_step(sdir, erate, inc, T)
        strain.append(np.dot(driver.strain_int[-1], sdir))
        stress.append(np.dot(driver.stress_int[-1], sdir))

  strain = np.array(strain)
  stress = np.array(stress)

  return {'strain': strain, 'stress': stress, 
      'energy_density': np.copy(driver.u),
      'plastic_work': np.copy(driver.p)}


def isochronous_curve(model, time, T = 300.0, emax = 0.05, srate = 1.0,
    ds = 10.0, max_cut = 4, nsteps = 250, history = None,
    check_dmg = False):
  """
    Generates an isochronous stress-strain curve at the given time and
    temperature.

    Parameters:
      model     material model
      time      relevant time

    Optional:
      T         temperature
      emax      maximum strain to attain on the curves
      srate     stress rate to use in the ramp part
      ds        stress increment along the curve
      max_cut   adaptive refinement
  """
  def strain(stress):
    res = creep(model, stress, srate, time, T = T, nsteps = nsteps,
        history = history, check_dmg = check_dmg)
    return res['tstrain'][-1], res['failed']

  strains = [0.0]
  stresses = [0.0]
  ncut = 0
  failed = False
  try:
    while strains[-1] < emax:
      target = stresses[-1] + ds
      try:
        enext, failed = strain(target)
        if failed:
          break
        stresses.append(target)
        strains.append(enext)
      except Exception:
        ncut += 1
        if ncut > max_cut:
          raise MaximumSubdivisions()
        ds /= 2
  except MaximumSubdivisions:
    # We were quite aggressive, so assume the curve goes flat
    stresses.append(stresses[-1])
    strains.append(emax)

  if failed:
    stresses.append(stresses[-1])
    strains.append(emax)

  # Now interpolate back the last strain point
  iff = inter.interp1d(strains, stresses)
  ls = iff(emax)
  strains[-1] = emax
  stresses[-1] = ls

  return {'strain': np.copy(strains), 'stress': np.copy(stresses)}

def offset_stress(e, s, eo = 0.2/100.0):
  """
    Helper function to generate yield stress from offset stress/strain data
    
    Parameters:
      e     strain data
      s     stress data

    Optional:
      eo    strain offset
  """
  iff = inter.interp1d(e, s)
  E = s[1] / e[1]
  
  eoff = opt.brentq(lambda e: iff(e) - E * (e - eo), 0.0,np.max(e))
  soff = iff(eoff)

  return soff

def gauss_points(n):
  """
    Helper for the below
  """
  xi, wi = leggauss(n)
  return 0.5*xi + 0.5, 0.5*wi

def bree(models, P, dT, T0 = 0.0, ncycles = 5, nsteps_up = 15,
    nsteps_cycle = 15, dt_load = 1.0, dt_cycle = 1.0, quadrature = 'midpoint',
    xi = None, Ai = None):
  """
    Solve a Bree problem using an approximate solution to Bree's integral
    equation.

    Parameters:
      models        vector of material models, one at each dl
      lengths       each dl
      P             primary stress
      dT            temperature difference

    Optional:
      T0            initial temperature
      ncycles       number of loading cycles to run
      nsteps_up     number of steps to load up in
      nsteps_cycle  number of steps per half cycle
      dt_load       time increment for loading
      dt_cycle      time increment for cycling
  """
  n = len(models)

  if ((xi is not None) and (Ai is None)) or ((Ai is not None) and (xi is None)):
    raise ValueError("If you specify one of the areas or positions you must"
        " specify both!")

  if xi is not None:
    xi = np.array(xi)
    Ai = np.array(Ai)
    if (np.max(xi) > 1) or (np.min(xi) < 0):
      raise ValueError("xi should go from zero to one!")
  else:
    if quadrature == 'midpoint':
      xi = (2.0 * (np.array(range(n)) + 1) - 1.0) / (2.0 * n)
      Ai = 1.0 / n
    elif quadrature == 'gauss':
      xi, Ai = gauss_points(n)
    else:
      raise ValueError("Unknown quadrature rule %s" % quadrature)

  # Do some housekeeping
  hist = [m.init_store() for m in models]
  mech_strain = np.zeros((n,))
  stresses = np.zeros((n,))
  temps = np.ones((n,)) * T0
  tn = 0.0

  # Update to the next increment
  def update(e, dTi, dt, t_n, hists_n, emechs_n, stresses_n, T_n):
    tempss = xi * dTi + T0
    t_np1 = t_n + dt

    stressess = np.zeros((n,))
    tangents = np.zeros((n,))
    strains = np.zeros((n,))
    work = np.zeros((n,))
    energy = np.zeros((n,))
    hists = []

    for i in range(n):
      m_i = models[i]
      emech_i = e - (tempss[i] - T0) * m_i.alpha(temps[i])
      strains[i] = emech_i
      s_i, h_i, A_i, u_i, p_i = m_i.update(
          emech_i, emechs_n[i], tempss[i], T_n[i], t_np1, t_n, stresses_n[i],
          hists_n[i], 0.0, 0.0)

      stressess[i] = s_i
      tangents[i] = A_i
      work[i] = p_i
      energy[i] = u_i
      hists.append(h_i)

    return stressess, tangents, hists, strains, temps, energy, work
  
  # Function to solve at each step
  def RJ(e, dTi, Pi, dt, t_n, hists_n, emechs_n, stresses_n, T_n):
    stresses, tangents, hists, strains, temps, energy, work = update(e, dTi, dt, t_n,
        hists_n, emechs_n, stresses_n, T_n)
    R = np.sum(stresses * Ai) - Pi
    A = np.sum(tangents * Ai)

    return R, A

  over_strain = [0.0]
  over_time = [0.0]
  over_energy = [0.0]
  over_work = [0.0]

  # Load up the model
  for Pi in np.linspace(0, P, nsteps_up+1)[1:]:
    sfn = lambda x: RJ(x, 0.0, Pi, dt_load, over_time[-1], hist, mech_strain,
        stresses, temps)
    e = scalar_newton(sfn, over_strain[-1])
    stresses, tangents, hist, mech_strain, temps, energy, work = update(
        e, 0.0, dt_load, over_time[-1], hist, mech_strain,
        stresses, temps) 
    over_strain.append(e)
    over_time.append(over_time[-1] + dt_load)
    over_energy.append(over_energy[-1] + np.sum(energy * Ai))
    over_work.append(over_work[-1] + np.sum(work * Ai))

  # Cycle the model
  for ci in range(ncycles):
    for dTi in np.linspace(0, dT, nsteps_cycle+1)[1:]:
      sfn = lambda x: RJ(x, dTi, P, dt_load, over_time[-1], hist, mech_strain,
          stresses, temps)
      e = scalar_newton(sfn, over_strain[-1])
      stresses, tangents, hist, mech_strain, temps, energy, work = update(
          e, dTi, dt_cycle, over_time[-1], hist, mech_strain,
          stresses, temps) 
      over_strain.append(e)
      over_time.append(over_time[-1] + dt_cycle)
      over_energy.append(over_energy[-1] + np.sum(energy * Ai))
      over_work.append(over_work[-1] + np.sum(work * Ai))

    for dTi in np.linspace(dT, 0, nsteps_cycle+1)[1:]:
      sfn = lambda x: RJ(x, dTi, P, dt_load, over_time[-1], hist, mech_strain,
          stresses, temps)
      e = scalar_newton(sfn, over_strain[-1])
      stresses, tangents, hist, mech_strain, temps, energy, work = update(
          e, dTi, dt_cycle, over_time[-1], hist, mech_strain,
          stresses, temps) 
      over_strain.append(e)
      over_time.append(over_time[-1] + dt_cycle)
      over_energy.append(over_energy[-1] + np.sum(energy * Ai))
      over_work.append(over_work[-1] + np.sum(work * Ai))

  return over_time, over_strain, over_energy, over_work


class MaximumIterations(RuntimeError):
  """
    An error to use if an iterative method exceeds the maximum allowed iterations.
  """
  def __init__(self):
    super(MaximumIterations, self).__init__("Exceeded the maximum allowed iterations!")

class MaximumSubdivisions(RuntimeError):
  """
    An error to use if adaptive substepping faisl.
  """
  def __init__(self):
    super(MaximumSubdivisions, self).__init__("Exceeded the maximum allowed step subdivisions!")


def newton(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-10, miter = 50,
    linesearch = 'none', bt_tau = 0.5, bt_c = 1.0e-4):
  """
    Manually-code newton-raphson so that I can output convergence info, if
    requested.

    Parameters:
      RJ        function return the residual + jacobian
      x0        initial guess

    Optional:
      verbose   verbose output
  """
  R, J = RJ(x0)
  nR = la.norm(R)
  nR0 = nR
  x = np.copy(x0)
  
  i = 0

  if verbose:
    print("Iter.\tnR\t\tnR/nR0\t\tcond\t\tlinesearch")
    print("%i\t%e\t%e\t" % (i, nR, nR / nR0))

  while (nR > rtol * nR0) and (nR > atol):
    a = la.solve(J, R)

    if linesearch == 'none':
      f = 1.0
    elif linesearch == 'backtracking':
      f = backtrack(RJ, R, J, x, -a, tau = bt_tau, c = bt_c, verbose = verbose)
    else:
      raise ValueError("Unknown linesearch type.")
    
    x -= (a * f)
    R, J = RJ(x)
    nR = la.norm(R)
    i += 1
    if verbose:
      print("%i\t%e\t%e\t%e\t%f" % (i, nR, nR / nR0,la.cond(J), f))
    if i > miter:
      if verbose:
        print("")
      raise MaximumIterations()
  
  if verbose:
    print("")

  return x

def backtrack(RJ, R0, J0, x, a, tau = 0.5, c = 1.0e-4, verbose = False):
  """
    Do a backtracking line search

    Parameters:
      RJ        residual/jacobian function
      R0        value of R, to save a feval
      J0        value of J, to save a feval
      x         point to start from
      a         direction

    Optional:
      tau       backtrack tau
      c         backtrack c
      verbose   verbose output
  """
  alpha = 1.0
  cv = la.norm(RJ(x + alpha * a)[0])
  while cv > la.norm(R0 + c * alpha * np.dot(J0, a)):
    alpha *= tau
    cv = la.norm(RJ(x + alpha * a)[0])

  return alpha

def quasi_newton(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-10, 
    miter = 20):
  """
    Manually-code quasi newton-raphson so that I can output convergence info, if
    requested.

    Parameters:
      RJ        function return the residual + jacobian
      x0        initial guess

    Optional:
      verbose   verbose output
  """
  R, J = RJ(x0)
  J_use = np.copy(J)
  nR = la.norm(R)
  nR0 = nR
  x = np.copy(x0)
  
  i = 0

  if verbose:
    print("Iter.\tnR\t\tnR/nR0\t\tcond")
    print("%i\t%e\t%e\t" % (i, nR, nR / nR0))

  while (nR > rtol * nR0) and (nR > atol):
    x -= la.solve(J_use, R)
    
    R, J = RJ(x)
    nR = la.norm(R)
    i += 1
    if verbose:
      print("%i\t%e\t%e\t%e" % (i, nR, nR / nR0,la.cond(J_use)))
    if i > miter:
      if verbose:
        print("")
      raise MaximumIterations()
  
  if verbose:
    print("")

  return x

def scalar_newton(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-10, miter = 20,
    quasi = None):
  """
    Manually-code newton-raphson so that I can output convergence info, if
    requested.

    Parameters:
      RJ        function return the residual + jacobian
      x0        initial guess

    Optional:
      verbose   verbose output
  """
  R, J = RJ(x0)
  nR = np.abs(R)
  nR0 = nR
  x = x0
  
  i = 0

  if verbose:
    print("Iter.\tnR\t\tnR/nR0")
    print("%i\t%e\t%e\t" % (i, nR, nR / nR0))

  while (nR > rtol * nR0) and (nR > atol):
    if quasi is not None:
      dx = -R/J
      x += quasi * dx
    else:
      x -= R/J
    R, J = RJ(x)
    nR = np.abs(R)
    i += 1
    if verbose:
      print("%i\t%e\t%e\t" % (i, nR, nR / nR0))
    
    if i > miter:
      if verbose:
        print("")
      raise MaximumIterations()
  
  if verbose:
    print("")

  return x
