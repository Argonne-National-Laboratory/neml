import numpy as np
import numpy.linalg as la

import scipy.interpolate as inter
import scipy.optimize as opt

class Driver(object):
  """
    Superclass of all drivers, basically just sets up history and reports
    results.
  """
  def __init__(self, model, verbose = False, rtol = 1.0e-6, atol = 1.0e-10,
      miter = 1000):
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

    self.stress_int = [np.zeros((6,))]

    self.stored_int = [self.model.init_store()]
    self.T_int = [0.0]
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

  @property
  def strain(self):
    return np.array(self.strain_int)

  def strain_step(self, e_np1, t_np1, T_np1):
    """
      Take a strain-controlled step

      Parameters:
        e_np1       next strain
        t_np1       next time
        T_np1       next temperature
    """
    s_np1, h_np1, A_np1, u_np1, p_np1 = self.model.update_sd(e_np1, 
        self.strain_int[-1],
        T_np1, self.T_int[-1], t_np1, self.t_int[-1], self.stress_int[-1],
        self.stored_int[-1], self.u_int[-1], self.p_int[-1])

    self.strain_int.append(np.copy(e_np1))
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
    def RJ(e):
      s, h, A, u, p = self.model.update_sd(e, self.strain_int[-1],
        T_np1, self.T_int[-1], t_np1, self.t_int[-1], 
        self.stress_int[-1],
        self.stored_int[-1], self.u_int[-1], self.p_int[-1])
      R = s - s_np1
      return R, A

    e_np1 = newton(RJ, self.strain_int[-1], verbose = self.verbose,
        rtol = self.rtol, atol = self.atol, miter = self.miter)

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
    def RJ(x):
      a = x[0]
      e_inc = x[1:]
      s, h, A, u, p = self.model.update_sd(self.strain_int[-1] + e_inc,
          self.strain_int[-1],
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

    x = newton(RJ, x0, verbose = self.verbose,
        rtol = self.rtol, atol = self.atol, miter = self.miter)
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
    s_np1 = self.stress_int[-1] + sdir / la.norm(sdir) * sinc
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
    oset = sorted(list(set(range(6)) - set([i])))
    def RJ(e_np1):
      s, h, A, u, p = self.model.update_sd(e_np1,
          self.strain_int[-1],
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
    e_np1 = newton(RJ, x0, verbose = self.verbose,
        rtol = self.rtol, atol = self.atol, miter = self.miter)

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
    super(Driver_sd_twobar, self).__init__(*args, **kwargs)

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

      em1 = e1 - np.array([1,1,1,0,0,0]) * a1 * (T1 - self.T1_0)
      em2 = e2 - np.array([1,1,1,0,0,0]) * a2 * (T2 - self.T2_0)
      
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
    try:
      x = newton(RJ, x0, verbose = self.verbose,
          rtol = self.rtol, atol = self.atol, miter = self.miter)
    except Exception:
      try:
        x = newton(RJ, x0, verbose = self.verbose,
            rtol = self.rtol, atol = self.atol, miter = self.miter,
            quasi = 0.2)  
      except Exception:
        x = newton(RJ, x0, verbose = False,
            rtol = self.rtol, atol = self.atol, miter = 10000,
            quasi = 1.0e-2)  


    e1 = x[:6]
    e2 = x[6:]

    self.strain1_int.append(e1)
    self.strain2_int.append(e2)

    a1 = self.model.alpha(T1)
    a2 = self.model.alpha(T2)

    em1 = e1 - np.array([1,1,1,0,0,0]) * a1 * (T1 - self.T1_0)
    em2 = e2 - np.array([1,1,1,0,0,0]) * a2 * (T2 - self.T2_0)

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

    estrain1 = self.model.elastic_strains(s1_np1, T1)
    estrain2 = self.model.elastic_strains(s2_np1, T2)

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
    dt = self.period / self.nsteps_cycle
    for toffset in np.linspace(0, self.period, self.nsteps_cycle+1)[1:]:
      self.take_step(self.P, self.T1_fn(toffset), self.T2_fn(toffset), dt,
          self.dstrain(toffset))

def twobar_test(model, A1, A2, T1, T2, period, P, load_time, ncycles,
    max_strain = None, min_strain = None, 
    nsteps_load = 50, nsteps_cycle = 201, verbose = False,
    rtol_classify = 1.0e-4, atol_classify = 1.0e-10,
    dstrain = lambda t: 0.0):
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
  """
  driver = Driver_sd_twobar(A1, A2, T1, T2, period, P, load_time, model,
      nsteps_load = nsteps_load, nsteps_cycle = nsteps_cycle, verbose = verbose,
      dstrain = dstrain)
  
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
        'classification': "collapse"
        }
  
  cycle_count = 0
  for i in range(ncycles):
    try:
      driver.one_cycle()
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
          'classification': "unstable"
          }
    ep1 = driver.strain1_plastic_int[-1][0]
    ep2 = driver.strain2_plastic_int[-1][0]
    max_s = np.max([np.abs(ep1), np.abs(ep2)])
    min_s = np.min([np.abs(ep1), np.abs(ep2)])
    cycle_count += 1
    if (max_strain is not None) and (max_s > max_strain):
      break
    if (min_strain is not None) and (min_s > min_strain):
      break

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
  
  return {
      'strain': driver.strain[nsteps_load:], 
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
      'strain2': driver.strain2[nsteps_load:,0]
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
    offset = 0.2/100.0):
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

    Results:
      strain    strain in direction
      stress    stress in direction
  """
  e_inc = emax / nsteps
  driver = Driver_sd(model, verbose = verbose)
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
  sYe = opt.brentq(lambda e: sfn(e) - tfn(e), 0.0, np.max(strain))
  sY = tfn(sYe)

  return {'strain': strain, 'stress': stress, 
      'energy_density': np.copy(driver.u),
      'plastic_work': np.copy(driver.p),
      'youngs': E, 'yield': sY}

def strain_cyclic(model, emax, R, erate, ncycles, T = 300.0, nsteps = 50,
    sdir = np.array([1,0,0,0,0,0]), hold_time = None, n_hold = 10,
    verbose = False):
  """
    HOLD IS BROKEN
    We need a generalized version of the strain_hold_step above to
    make it work, but I'm tired right now.

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
  driver = Driver_sd(model, verbose = verbose)
  emin = emax * R
  if hold_time:
    raise NotImplementedError("Whoops, this is broken")
    if np.isscalar(hold_time):
      hold_time = [hold_time, hold_time]

  # Setup results
  strain = [0.0]
  stress = [0.0]
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
  for i in range(nsteps):
    if i == 0:
      einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T)
    else:
      einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T, einc_guess = einc,
          ainc_guess = ainc)
    strain.append(np.dot(driver.strain_int[-1], sdir))
    stress.append(np.dot(driver.stress_int[-1], sdir))
  
  # Begin cycling
  for s in range(ncycles):
    if verbose:
      print("Cycle %i" % s)

    si = len(driver.strain_int)
    e_inc = np.abs(emin - emax) / nsteps
    for i in range(nsteps):
      if i == 0:
        einc, ainc = driver.erate_einc_step(-sdir, erate, e_inc, T, 
            einc_guess = -einc, ainc_guess = -ainc)
      else:
        einc, ainc = driver.erate_einc_step(-sdir, erate, e_inc, T, 
            einc_guess = einc, ainc_guess = ainc)

      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

    e_inc = np.abs(emax - emin) / nsteps
    for i in range(nsteps):
      if i == 0:
        einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T, 
            einc_guess = -einc, ainc_guess = -ainc)
      else:
        einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T,
            einc_guess = einc, ainc_guess = ainc)
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

    # Calculate
    cycles.append(s)
    smax.append(max(stress[si:]))
    smin.append(min(stress[si:]))
    smean.append((smax[-1]+smin[-1])/2)

    ecycle.append(driver.u_int[-1])
    pcycle.append(driver.p_int[-1])

  # Setup and return
  return {"strain": np.array(strain), "stress": np.array(stress),
      "cycles": np.array(cycles, dtype = int), "max": np.array(smax),
      "min": np.array(smin), "mean": np.array(smean),
      "energy_density": np.array(ecycle), "plastic_work": np.array(pcycle),
      "history": driver.stored_int[-1]}

def stress_cyclic(model, smax, R, srate, ncycles, T = 300.0, nsteps = 50,
    sdir = np.array([1,0,0,0,0,0]), hold_time = None, n_hold = 10,
    verbose = False):
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
  driver = Driver_sd(model, verbose = verbose)
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
  driver = Driver_sd(model, verbose = verbose)
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


def creep(model, smax, srate, hold, T = 300.0, nsteps = 750,
    nsteps_up = 150, sdir = np.array([1,0,0,0,0,0]), verbose = False,
    logspace = False, history = None):
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
  """
  # Setup
  driver = Driver_sd(model, verbose = verbose)
  if history is not None:
    driver.stored_int[0] = history
  time = [0]
  strain = [0]
  stress = [0]

  # Ramp up
  sinc = smax / nsteps_up
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
   
  for t in ts:
    # You can exceed the creep life of the sample doing this...
    # Need to allow a non-convergent result to break
    try:
      driver.stress_step(driver.stress_int[-1], t, T)
    except:
      break
    time.append(t)
    strain.append(np.dot(driver.strain_int[-1],sdir))
    stress.append(np.dot(driver.stress_int[-1],sdir))

  time = np.array(time)
  strain = np.array(strain)
  stress = np.array(stress)
  rrate = np.diff(strain[ri:]) / np.diff(time[ri:])
  rstrain = strain[ri:]
  
  return {'time': np.copy(time), 'strain': np.copy(strain), 
      'stress': np.copy(stress), 'rtime': np.copy(time[ri:-1] - time[ri]),
      'rrate': np.copy(rrate), 'rstrain': np.copy(rstrain[:-1])}

def rate_jump_test(model, erates, T = 300.0, e_per = 0.01, nsteps_per = 100, 
    sdir = np.array([1,0,0,0,0,0]), verbose = False, history = None):
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

    Results:
      strain    strain in direction
      stress    stress in direction
  """
  e_inc = e_per / nsteps_per
  driver = Driver_sd(model, verbose = verbose)
  if history is not None:
    driver.stored_int[0] = history
  strain = [0.0]
  stress = [0.0]
  
  for erate in erates:
    for i in range(nsteps_per):
      if i == 0:
        einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T)
      else:
        einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T, 
            einc_guess = einc, ainc_guess = ainc)
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

  strain = np.array(strain)
  stress = np.array(stress)

  return {'strain': strain, 'stress': stress, 
      'energy_density': np.copy(driver.u),
      'plastic_work': np.copy(driver.p)}


def isochronous_curve(model, time, T = 300.0, emax = 0.05, srate = 1.0e-2,
    ds = 10.0, max_cut = 10, nsteps = 250, history = None):
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
        history = history)
    return res['rstrain'][-1]

  strains = [0.0]
  stresses = [0.0]
  ncut = 0
  try:
    while strains[-1] < emax:
      #print(strains[-1])
      target = stresses[-1] + ds
      try:
        enext = strain(target)
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
  E = s[3] / e[3]
  
  eoff = opt.brentq(lambda e: iff(e) - E * (e - eo), 0.0,np.max(e))
  soff = iff(eoff)

  return soff

def bree(models, lengths, P, dT, T0 = 0.0, ncycles = 5, nsteps_up = 15,
    nsteps_cycle = 15, dt_load = 1.0, dt_cycle = 1.0):
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
  # Do some housekeeping
  lengths = np.array(lengths)
  n = len(models)
  hist = [m.init_store() for m in models]
  mech_strain = np.zeros((n,))
  stresses = np.zeros((n,))
  temps = np.ones((n,)) * T0
  tn = 0.0

  # Centers for each segment
  ends = np.cumsum(lengths)
  xs = ends - np.array(lengths) / 2.0
  L = np.sum(lengths)

  # Update to the next increment
  def update(e, dTi, dt, t_n, hists_n, emechs_n, stresses_n, T_n):
    tempss = xs / L * dTi + T0
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
    R = np.sum(stresses * lengths) - Pi * L
    A = np.sum(tangents * lengths)

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
    over_energy.append(over_energy[-1] + np.sum(energy * lengths))
    over_work.append(over_work[-1] + np.sum(work * lengths))

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
      over_energy.append(over_energy[-1] + np.sum(energy * lengths))
      over_work.append(over_work[-1] + np.sum(work * lengths))

    for dTi in np.linspace(dT, 0, nsteps_cycle+1)[1:]:
      sfn = lambda x: RJ(x, dTi, P, dt_load, over_time[-1], hist, mech_strain,
          stresses, temps)
      e = scalar_newton(sfn, over_strain[-1])
      stresses, tangents, hist, mech_strain, temps, energy, work = update(
          e, dTi, dt_cycle, over_time[-1], hist, mech_strain,
          stresses, temps) 
      over_strain.append(e)
      over_time.append(over_time[-1] + dt_cycle)
      over_energy.append(over_energy[-1] + np.sum(energy * lengths))
      over_work.append(over_work[-1] + np.sum(work * lengths))

  return over_time, over_strain, over_energy, over_work


class MaximumIterations(RuntimeError):
  def __init__(self):
    super(MaximumIterations, self).__init__("Exceeded the maximum allowed iterations!")

class MaximumSubdivisions(RuntimeError):
  def __init__(self):
    super(MaximumSubdivisions, self).__init__("Exceeded the maximum allowed step subdivisions!")


def newton(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-10, miter = 20,
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
  nR = la.norm(R)
  nR0 = nR
  x = np.copy(x0)
  
  i = 0

  if verbose:
    print("Iter.\tnR\t\tnR/nR0")
    print("%i\t%e\t%e\t" % (i, nR, nR / nR0))

  while (nR > rtol * nR0) and (nR > atol):
    if quasi is not None:
      dx = -la.solve(J, R)
      x += quasi * dx
    else:
      x -= la.solve(J, R)
    R, J = RJ(x)
    nR = la.norm(R)
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
