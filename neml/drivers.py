import numpy as np
import numpy.linalg as la

import scipy.interpolate as inter
import scipy.optimize as opt
from numpy.polynomial.legendre import leggauss

import numpy.random as ra

from neml.nlsolvers import MaximumIterations, MaximumSubdivisions, newton, scalar_newton

class Driver(object):
  """
    Superclass of all drivers, basically just sets up history and reports
    results.
  """
  def __init__(self, model, verbose = False, rtol = 1.0e-6, atol = 1.0e-10,
      miter = 25, T_init = 0.0, no_thermal_strain = False):
    """
      Parameters:
        model:       material model to play with
        verbose:     verbose output
        rtol:        relative tolerance, where needed
        atol:        absolute tolerance, where needed
        miter:       maximum iterations, where needed
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
    """
      Parameters:
        model:       material model to play with
        verbose:     verbose output
        rtol:        relative tolerance, where needed
        atol:        absolute tolerance, where needed
        miter:       maximum iterations, where needed
    """
    super(Driver_sd, self).__init__(*args, **kwargs)
    self.strain_int = [np.zeros((6,))]
    self.thermal_strain_int = [np.zeros((6,))]
    self.mechanical_strain_int = [np.zeros((6,))]

  def solve_try(self, RJ, x0, extra = []):
    """
      Try several different nonlinear solvers in the hope that at least
      one will converge

      Parameters:
        RJ:      function that returns the residual equations and associated
                 Jacobian
        x0:      initial guess
        extra:   list of extra solver functions of the type below
    """
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

    solvers = [s1,s3]
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
    """
      Move the thermal strains to the next step

      Parameters:
        T_np1:      next temperature
    """
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
        e_np1:       next strain
        t_np1:       next time
        T_np1:       next temperature
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
        s_np1:       next stress
        t_np1:       next time
        T_np1:       next temperature
    """
    enext = self.update_thermal_strain(T_np1)
    def RJ(e):
      s, h, A, u, p = self.model.update_sd(e - enext, self.mechanical_strain_int[-1],
        T_np1, self.T_int[-1], t_np1, self.t_int[-1],
        self.stress_int[-1],
        self.stored_int[-1], self.u_int[-1], self.p_int[-1])
      R = s - s_np1
      return R, A

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
        sdir:        stress direction
        erate:       strain rate (in the direction)
        t_np1:       next time
        T_np1:       next temperature
        einc_guess:  a guess at the strain increment
        ainc_guess:  a guess at the stress increment
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

    x = self.solve_try(RJ, x0)
    e_np1 = self.strain_int[-1] + x[1:]

    self.strain_step(e_np1, t_np1, T_np1)

    return x[1:], x[0]

  def erate_einc_step(self, sdir, erate, einc, T_np1, **kwargs):
    """
      Similar to erate_step but specify the strain increment instead of the
      time increment.

      Parameters:
        sdir:        stress direction
        erate:       strain rate, in stress direction
        einc:        strain increment, in stress direction
        T_np1:       temperature at next time step
    """
    dt = einc / erate
    return self.erate_step(sdir, erate, self.t_int[-1] + dt, T_np1, **kwargs)

  def srate_sinc_step(self, sdir, srate, sinc, T_np1):
    """
      Similar to rate_step but specify the stress increment instead of the
      time increment.

      Parameters:
        sdir:        stress direction
        srate:       stress rate
        sinc:        stress increment
        T_np1:       temperature at next time step

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

  def strain_hold_step(self, i, t_np1, T_np1, q = 1.0, E = -1.0):
    """
      A special, mixed step which holds the strain in index i constant
      while holding the stress in the other directions to their previous
      values

      Parameters:
        i:           index to hold
        t_np1:       next time
        T_np1:       next temperature
        q:           follow up factor
        E:           Young's modulus to use -- must redo interface at some point
    """
    if not np.isclose(q, 1.0) and np.isclose(E, -1.0):
      raise ValueError("You must supply the Youngs modulus")

    enext = self.update_thermal_strain(T_np1)
    oset = sorted(list(set(range(6)) - set([i])))
    def RJ(e_np1):
      s, h, A, u, p = self.model.update_sd(e_np1 - enext,
          self.mechanical_strain_int[-1],
          T_np1, self.T_int[-1], t_np1, self.t_int[-1], self.stress_int[-1],
          self.stored_int[-1], self.u_int[-1], self.p_int[-1])

      R = np.zeros((6,))
      R[0] = (e_np1[i] - self.strain_int[-1][i]
          ) + (s[i] - self.stress_int[-1][i]) / E * (q - 1)
      R[1:] = s[oset] - self.stress_int[-1][oset]

      J = np.zeros((6,6))
      J[0,0] = 1.0
      J[0,:] += A[i,:] / E * (q - 1)
      J[1:,:] = A[oset,:][:]

      return R, J

    x0 = np.copy(self.strain_int[-1])

    e_np1 = self.solve_try(RJ, x0)

    self.strain_step(e_np1, t_np1, T_np1)

def uniaxial_test(model, erate, T = 300.0, emax = 0.05, nsteps = 250,
    sdir = np.array([1,0,0,0,0,0]), verbose = False,
    offset = 0.2/100.0, history = None, tdir = np.array([0,1,0,0,0,0]),
    rtol = 1e-6, atol = 1e-10, miter = 25,
    full_results = False):
  """
    Make a uniaxial stress/strain curve

    Parameters:
      model:            material model
      erate:            strain rate

    Keyword Args:
      T:                temperature, default 300.0
      emax:             maximum strain, default 5%
      nsteps:           number of steps to use, default 250
      sdir:             stress direction, default tension in x
      verbose:          whether to be verbose
      offset:           used to calculate yield stress
      history:          initial model history
      tdir:             transverse direction for Poisson's ratio

    Returns:
      dict:             results dictionary containing...


    **Results in dictionary:**
      ================= ============================================
      Name              Description
      ================= ============================================
      strain            strain in direction
      stress            stress in direction
      energy_density    strain energy density
      plastic_work      plastic dissipation
      youngs            young's modulus of initial curve
      yield             yield stress implied by curve
      poissons          poisson's ratio implied by non-axial strains
      ================= ============================================
  """
  e_inc = emax / nsteps
  driver = Driver_sd(model, verbose = verbose, T_init = T, rtol = rtol,
      atol = atol, miter = miter)
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
  nu = -np.dot(driver.strain_int[1], tdir) / np.dot(
      driver.strain_int[1], sdir)
  sfn = inter.interp1d(np.abs(strain), np.abs(stress))
  tfn = lambda e: E * (e - offset)
  try:
    sYe = opt.brentq(lambda e: sfn(e) - tfn(e), 0.0, np.max(strain))
    sY = tfn(sYe)
  except Exception:
    sY = np.inf
  
  if full_results:
    return {'time': driver.t_int,
        'strain': driver.mechanical_strain_int,
        'temperature': driver.T_int,
        'stress': driver.stress_int }
  else:
    return {'strain': strain, 'stress': stress,
        'energy_density': np.copy(driver.u),
        'plastic_work': np.copy(driver.p),
        'youngs': E, 'yield': sY, 'poissons': nu,
        'history': driver.stored_int[-1]}

def strain_cyclic(model, emax, R, erate, ncycles, T = 300.0, nsteps = 50,
    sdir = np.array([1,0,0,0,0,0]), hold_time = None, n_hold = 25,
    verbose = False, check_dmg = False, dtol = 0.75):
  """
    Strain controlled cyclic test.

    Parameters:
      emax:         maximum strain
      R:            R = emin / emax
      erate:        strain rate to go at
      ncycles:      number of cycles
      T:            temperature, default 300

    Keyword Args:
      nsteps:       number of steps per half cycle
      sdir:         stress direction, defaults to x and tension first
      hold_time:    if None don't hold, if scalar then hold symmetrically top/bot
                    if an array specify different hold times for first direction
                    (default tension) and second direction
      n_hold:       number of steps to hold over
      verbose:      whether to be verbose
      check_dmg:    check to see if material damage exceeds dtol, stop the
                    simulation when that happens
      dtol:         damage to stop at

    Returns:
      dict:         results dictionary containing...

    **Results in dictionary:**
      ============= ========================
      Name          Description
      ============= ========================
      strain:       strain in direction
      stress:       stress in direction
      cycles:       list of cycle numbers
      max:          maximum stress per cycle
      min:          minimum stress per cycle
      mean:         mean stress per cycle
      ============= ========================
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

def strain_cyclic_extrapolated(model, emax, R, erate, ncycles, T = 300.0, nsteps = 50,
    sdir = np.array([1,0,0,0,0,0]), hold_time = None, n_hold = 25,
    verbose = False, check_dmg = False, dtol = 0.75, min_cycle=3, unit_extrapolate = 10,
    jump_delta_N=10, allowable_jump_stress=5.0):
  """
    Strain controlled cyclic test extrapolation.

    Extra Keyword Args:
      min_cycle                minimum cycles to start the extrapolation process
      unit_extrapolate         number of cycles to perform single cycle extrapolation
      jump_delta_N               number of cycles to jump
      allowable_jump_stress    extrapolate when stress jump is within this limit

    Returns:
      dict:         results dictionary containing...

    **Results in dictionary:**
      ============= ========================
      Name          Description
      ============= ========================
      cycles:       list of cycle numbers
      max:          maximum stress per cycle
      min:          minimum stress per cycle
      ============= ========================
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
  s = 0

# steps in one cycle
  if (hold_time[0] > 0) and (hold_time[1] == 0):
      steps = 2*nsteps + n_hold
  elif (hold_time[1] > 0) and (hold_time[0] == 0):
      steps = 2*nsteps + n_hold
  elif (hold_time[0] > 0) and (hold_time[1] > 0):
      steps = 2*nsteps + 2*n_hold
  else:
      steps = 2*nsteps

  extrapolate = False
  while s < ncycles:
    if verbose:
      print("Cycle %i" % s)

    if check_dmg:
      if driver.stored_int[-1][0] > dtol:
        print("Damage check exceeded")
        break

    if (s >= min_cycle) and (extrapolate == True):        # No extrapolation before min_cycle
        if (s <= unit_extrapolate):                       # single cycle jump for first unit_extrapolate cycles
            delta_N = 1
        else:
            delta_N = jump_delta_N                        # specified cycles to jump
        n = len(driver.stored_int)
        # extrapolating history
        pos_hist_last_last = driver.stored_int[n - 1 - steps]
        pos_hist_last      = driver.stored_int[n-1]
        dN_1 = cycles[-1] - cycles[-2]
        pos_extrapolated_history = pos_hist_last + (pos_hist_last - pos_hist_last_last)*delta_N/dN_1
        # extrapolating smax
        smax_last_last = smax[-2]
        smax_last = smax[-1]
        extrapolated_smax = smax_last + (smax_last - smax_last_last)*delta_N/dN_1
        # extrapolating smax
        smin_last_last = smin[-2]
        smin_last = smin[-1]
        extrapolated_smin = smin_last + (smin_last - smin_last_last)*delta_N/dN_1
        # criteria for extrapolation
        pos_stress_last_last = driver.stress_int[n - 1 - 2*steps]
        pos_stress_last      = driver.stress_int[n-1]
        pos_extrapolated_stress = pos_stress_last + (pos_stress_last - pos_stress_last_last)*delta_N/dN_1
        stress_jump = pos_extrapolated_stress[0] - pos_stress_last[0]
        if np.fabs(stress_jump) <= allowable_jump_stress:
            s = s + delta_N
            if s > ncycles:
                break
            driver.stored_int.append(pos_extrapolated_history)
            driver.stress_int.append(pos_extrapolated_stress)
            smax.append(extrapolated_smax)
            smin.append(extrapolated_smin)
            cycles.append(s)
            extrapolate = False
        else:
            extrapolate = False

    else:
        try:
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
          s += 1
          cycles.append(s)
          smax.append(max(stress[si:]))
          smin.append(min(stress[si:]))
          smean.append((smax[-1]+smin[-1])/2)
          ecycle.append(driver.u_int[-1])
          pcycle.append(driver.p_int[-1])
          extrapolate = True
        except Exception as e:
          break

   # Setup and return
  return {"cycles": np.array(cycles, dtype = int), "max": np.array(smax),
       "min": np.array(smin),"time": np.array(time)}


def strain_cyclic_followup(model, emax, R, erate, ncycles,
    q = 1.0, T = 300.0, nsteps = 50,
    sind = 0, hold_time = None, n_hold = 25,
    verbose = False, check_dmg = False, dtol = 0.75,
    logspace = False):
  """
    Strain controlled cyclic test with follow up.

    This is a "fallback" to the old version that does things by index
    so that I can use the index-based hold routine with follow up

    Parameters:
      emax:      maximum strain
      R:         R = emin / emax
      erate:     strain rate to go at
      ncycles:   number of cycles

    Keyword Args:
      q:         follow up factor
      T:         temperature, default 300
      nsteps:    number of steps per half cycle
      sind:      index to pull on
      hold_time: if None don't hold, if scalar then hold symmetrically top/bot
                 if an array specify different hold times for first direction
                 (default tension) and second direction
      n_hold:    number of steps to hold over
      verbose:   whether to be verbose
      check_dmg: check to see if damage exceeds a threshold
      dtol:      damage threshold
      logspace:  logspace the hold time steps (instead of linspace)

    Returns:
      dict:      dictionary of results...


    **Results in dictionary**
      ========= ========================
      Name      Description
      ========= ========================
      strain    strain in direction
      stress    stress in direction
      cycles    list of cycle numbers
      max       maximum stress per cycle
      min       minimum stress per cycle
      mean      mean stress per cycle
      ========= ========================
  """
  # Setup
  sdir = np.zeros((6,))
  sdir[sind] = 1.0

  res = uniaxial_test(model, erate, T = T, emax = 1.0e-4, nsteps = 2)
  E = res['youngs']

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
        if logspace:
          dts = np.diff(np.logspace(0, np.log10(hold_time[0]), n_hold+1))
        else:
          dts = np.diff(np.linspace(0,hold_time[0],n_hold+1))
        #dt = hold_time[0] / n_hold
        for i, dt in enumerate(dts):
          driver.strain_hold_step(sind, time[-1] + dt, T,
              q = q, E = E)
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
              einc_guess = np.zeros((6,)), ainc_guess = -1)
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
        if logspace:
          dts = np.diff(np.logspace(0, np.log10(hold_time[1]), n_hold+1))
        else:
          dts = np.diff(np.linspace(0,hold_time[1],n_hold+1))
        for i, dt in enumerate(dts):
          driver.strain_hold_step(sind, time[-1] + dt, T,
              q = q, E = E)
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
              einc_guess = np.zeros((6,)), ainc_guess = 1.0)
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
      emax:         maximum stress
      R:            R = smin / smax
      erate:        strain rate to go at
      ncycles:      number of cycles

    Keyword Args:
      T:            temperature, default 300
      nsteps:       number of steps per half cycle
      sdir:         stress direction, defaults to x and tension first
      hold_time:    if None don't hold, if scalar then hold symmetrically top/bot
                    if an array specify different hold times for first direction
                    (default tension) and second direction
      n_hold:       number of steps to hold over
      verbose:      whether to be verbose
      etol:         stop if the strain increment exceeds this threshold

    Returns:
      dict:         dictionary of results


    **Results in dictionary:**
      ============= ========================
      Name          Description
      ============= ========================
      strain        strain in direction
      stress        stress in direction
      cycles        list of cycle numbers
      max           maximum strain per cycle
      min           minimum strain per cycle
      mean          mean strain per cycle
      ============= ========================
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
    if hold_time and (hold_time[0] > 0.0):
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

def stress_relaxation(model, emax, erate, hold, T = 300.0, nsteps = 250,
    nsteps_up = 50, index = 0, tc = 1.0,
    verbose = False, logspace = False, q = 1.0):
  """
    Simulate a stress relaxation test.

    Parameters:
      model:         material model
      emax :         maximum strain to attain
      erate:         strain rate to take getting there
      hold:          hold time

    Keyword Args:
      T:             temperature
      nsteps:        number of steps to relax over
      nsteps_up:     number of steps to take getting up to stress
      index:         direction to pull in, default x tension
      tc:            1.0 for tension -1.0 for compression
      verbose:       whether to be verbose
      logspace:      log space the relaxation timesteps
      q:             follow up factor

    Results:
      dict:          dictionary of results

    **Results in dictionary:**
      ============== ======================
      Name           Description
      ============== ======================
      time           time
      strain         strain
      stress         stress
      rtime          relaxation time
      rrate          stress relaxation rate
      ============== ======================
  """
  # Setup
  driver = Driver_sd(model, verbose = verbose, T_init = T)
  time = [0]
  strain = [0]
  stress = [0]

  res = uniaxial_test(model, erate, T = T, emax = 1.0e-4, nsteps = 2)
  E = res['youngs']

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
  if logspace:
    ts = np.logspace(0, np.log10(hold), num = nsteps+1)
    dts = np.diff(ts)
  else:
    dt = hold / nsteps
    dts = [dt] * nsteps
  for i, dt in enumerate(dts):
    driver.strain_hold_step(index, driver.t_int[-1] + dt, T,
        q = q, E = E)
    time.append(driver.t_int[-1])
    strain.append(np.dot(driver.strain_int[-1],sdir))
    stress.append(np.dot(driver.stress_int[-1],sdir))

  time = np.array(time)
  strain = np.array(strain)
  stress = np.array(stress)
  rrate = -np.diff(stress[ri:]) / np.diff(time[ri:])

  return {'time': np.copy(time), 'strain': np.copy(strain),
      'stress': np.copy(stress), 'rtime': np.copy(time[ri:-1] - time[ri]),
      'rrate': np.copy(rrate), 'rstress': np.copy(stress[ri:-1]),
      'rstrain': np.copy(strain[ri:-1])}

def creep(model, smax, srate, hold, T = 300.0, nsteps = 250,
    nsteps_up = 150, sdir = np.array([1,0,0,0,0,0]), verbose = False,
    logspace = False, history = None, elimit = 1.0, check_dmg = False,
    dtol = 0.75):
  """
    Simulate a creep test

    Parameters:
      model:         material model
      smax:          stress to attain
      srate:         stress rate to take getting there
      hold:          total hold time

    Keyword Args:
      T:             temperature
      nsteps:        number of steps over relaxation period
      nsteps_up:     number of steps to get to stress value
      sdir:          stress direction, defaults to x-tension
      verbose:       whether to be verbose
      logspace:      if true logspace the time steps
      history:       use damaged material
      check_dmg:     check damage as a break condition
      dtol:          damage to define failure at

    Returns:
      dict:          results dictionary
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
      'tstrain': np.copy(strain[ri:-1]),
      'history': np.array(driver.stored_int), 'failed': failed}

def thermomechanical_strain_raw(model, time, temperature, strain,
    sdir = np.array([1,0,0,0,0,0.0]), verbose = False, substep = 1):
  """
    Directly drive a model using the output of a strain controlled
    thermomechanical test

    Parameters:
      model:         material model
      time:          list of times
      temperature:   list of temperatures
      strain:        list of strains

    Keyword Args:
      sdir:          direction of stress
      verbose:       print information
      substep:       take substep steps per data point

    Returns:
      dict:          results dictionary
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
      model:         material model
      erate:         list of strain rates

    Keyword Args:
      T:             temperature, default 300.0
      e_per:         how much straining to do for each rate
      nsteps_per:    number of steps per strain rate
      sdir:          stress direction, default tension in x
      verbose:       whether to be verbose
      history:       prior model history
      strains:       manual definition of jumps

    Returns:
      dict:          results dictionary
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
    check_dmg = False, dtol = 0.75):
  """
    Generates an isochronous stress-strain curve at the given time and
    temperature.

    Parameters:
      model:     material model
      time:      relevant time

    Keyword Args:
      T:         temperature
      emax:      maximum strain to attain on the curves
      srate:     stress rate to use in the ramp part
      ds:        stress increment along the curve
      max_cut:   adaptive refinement
      nsteps:    number of creep steps
      history:   start with a non-zero initial history
      check_dmg: stop if damage exceeds a threshold
      dtol:      damage threshold
  """
  def strain(stress):
    res = creep(model, stress, srate, time, T = T, nsteps = nsteps,
        history = history, check_dmg = check_dmg, dtol = dtol)
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
      e:     strain data
      s:     stress data
      eo:    strain offset
  """
  iff = inter.interp1d(e, s)
  E = s[1] / e[1]

  eoff = opt.brentq(lambda e: iff(e) - E * (e - eo), 0.0,np.max(e))
  soff = iff(eoff)

  return soff

def classify(ua, ub, pa, pb, e1a, e1b, e2a, e2b, rtol = 1.0e-4, atol = 1.0e-10):
  """
    Classify a model as elastic, elastic shakedown, plastic shakedown,
    or ratcheting.

    Parameters:
      ua:    cycle a internal energy
      ub:    cycle b internal energy
      pa:    cycle a plastic dissipation
      pb:    cycle b plastic dissipation
      ea:    cycle a strain
      eb:    cycle b strain
      rtol:  relative tolerance
      atol:  absolute tolerance
  """
  if np.abs(pb) < atol:
    return 'elastic'
  elif np.abs(ub-ua) < rtol * np.abs(ub):
    return 'elastic shakedown'
  elif (np.abs(e1b-e1a) < rtol * np.abs(e1b)) and (np.abs(e2b-e2a) < rtol * np.abs(e2b)):
    return 'plastic shakedown'
  else:
    return 'ratcheting'

def def_grad_driver(model, F, tmax, nsteps, T = 300.0):
  """
    Basic large deformation driver pushing a model through the deformation
    gradient as a function of time

    Parameters:
      model     model to use
      F         deformation gradient as a function of time
      tmax      maximum time
      nsteps    number of load steps to take

    Optional:
      T         temperature
  """
  time = [0.0]
  stress = [np.zeros((6,))]
  deform = [F(0.0)]
  D_hist = [np.zeros((6,))]
  W_hist = [np.zeros((3,))]
  history = [model.init_store()]
  energy = [0.0]
  dissipation = [0.0]

  t = 0.0
  dt = tmax / nsteps

  for i in range(nsteps):
    t += dt
    F_np1 = F(t)
    L = np.dot((F_np1 - deform[-1]), la.inv(F_np1))
    dD = sym(0.5*(L + L.T))
    dW = skew(0.5*(L - L.T))

    s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = model.update_ld_inc(
        D_hist[-1] + dD, D_hist[-1],
        W_hist[-1] + dW, W_hist[-1],
        T, T, t, time[-1], stress[-1],
        history[-1], energy[-1], dissipation[-1])

    time.append(t)
    stress.append(s_np1)
    deform.append(F_np1)
    D_hist.append(D_hist[-1] + dD)
    W_hist.append(W_hist[-1] + dW)
    history.append(h_np1)
    energy.append(u_np1)
    dissipation.append(p_np1)

  return {'time': time, 'stress': stress, 'F': deform, 'D': D_hist,
      'W': W_hist, 'history': history, 'energy': energy,
      'dissipation': dissipation}

def sym(A):
  """
    Take a symmetric matrix to the Mandel convention vector.
  """
  return np.array([A[0,0], A[1,1], A[2,2], np.sqrt(2)*A[1,2],
    np.sqrt(2)*A[0,2], np.sqrt(2)*A[0,1]])

def skew(A):
  """
    Take a skew matrix to my vector convention.
  """
  return np.array([-A[1,2], A[0,2], -A[0,1]])
