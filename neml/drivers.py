import numpy as np
import numpy.linalg as la

import scipy.interpolate as inter

class Driver(object):
  """
    Superclass of all drivers, basically just sets up history and reports
    results.
  """
  def __init__(self, model, verbose = False, rtol = 1.0e-6, atol = 1.0e-10,
      miter = 20):
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


def uniaxial_test(model, erate, T = 300.0, emax = 0.05, nsteps = 250, 
    sdir = np.array([1,0,0,0,0,0]), verbose = False):
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

  return {'strain': strain, 'stress': stress, 
      'energy_density': np.copy(driver.u),
      'plastic_work': np.copy(driver.p)}

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

    si = len(driver.strain)
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
      "energy_density": np.array(ecycle), "plastic_work": np.array(pcycle)}

def stress_cyclic(model, smax, R, srate, ncycles, T = 300.0, nsteps = 50,
    sdir = np.array([1,0,0,0,0,0]), hold_time = None, n_hold = 10,
    verbose = False):
  """
    Strain controlled cyclic test.

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
    if verbose:
      print("Cycle %i" % s)
    si = len(driver.strain)
    # Hold, if requested
    if hold_time and (hold_time[0] > 0.0):
      ht = hold_time[0]
      dt = ht / n_hold
      for i in range(n_hold):
        driver.stress_step(driver.stress_int[-1], driver.t_int[-1] + dt, T)
        strain.append(np.dot(driver.strain_int[-1], sdir))
        stress.append(np.dot(driver.stress_int[-1], sdir))

    s_inc = (smin - smax) / nsteps
    for i in range(nsteps):
      driver.srate_sinc_step(sdir, srate, s_inc, T)
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

    # Hold, if requested
    if hold_time and (hold_time[1] > 0.0):
      ht = hold_time[1]
      dt = ht / n_hold
      for i in range(n_hold):
        driver.stress_step(driver.stress_int[-1], driver.t_int[-1] + dt, T)
        strain.append(np.dot(driver.strain_int[-1], sdir))
        stress.append(np.dot(driver.stress_int[-1], sdir))

    s_inc = (smax - smin) / nsteps
    for i in range(nsteps):
      driver.srate_sinc_step(sdir, srate, s_inc, T)
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

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
      "energy_density": np.array(ecycle), "plastic_work": np.array(pcycle)}

def stress_relaxation(model, emax, erate, hold, T = 300.0, nsteps = 500,
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


def creep(model, smax, srate, hold, T = 300.0, nsteps = 250,
    nsteps_up = 150, sdir = np.array([1,0,0,0,0,0]), verbose = False,
    logspace = False):
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
  """
  # Setup
  driver = Driver_sd(model, verbose = verbose)
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
    driver.stress_step(driver.stress_int[-1], t, T)
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
      'rrate': np.copy(rrate), 'rstrain': np.copy(rstrain)}

def isochronous_curve(model, time, T = 300.0, emax = 0.05, srate = 1.0e-2,
    ds = 10.0, max_cut = 10, nsteps = 250):
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
    res = creep(model, stress, srate, time, T = T, nsteps = nsteps)
    return res['rstrain'][-1]

  strains = [0.0]
  stresses = [0.0]
  ncut = 0
  try:
    while strains[-1] < emax:
      print(strains[-1])
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

class MaximumIterations(RuntimeError):
  def __init__(self):
    super(MaximumIterations, self).__init__("Exceeded the maximum allowed iterations!")

class MaximumSubdivisions(RuntimeError):
  def __init__(self):
    super(MaximumSubdivisions, self).__init__("Exceeded the maximum allowed step subdivisions!")


def newton(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-10, miter = 20):
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
