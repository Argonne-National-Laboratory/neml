import numpy as np
import numpy.linalg as la

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
    s_np1, h_np1, A_np1 = self.model.update_sd(e_np1, self.strain_int[-1],
        T_np1, self.T_int[-1], t_np1, self.t_int[-1], self.stress_int[-1],
        self.stored_int[-1])

    self.strain_int.append(np.copy(e_np1))
    self.stress_int.append(np.copy(s_np1))
    self.stored_int.append(np.copy(h_np1))
    self.T_int.append(T_np1)
    self.t_int.append(t_np1)

  def stress_step(self, s_np1, t_np1, T_np1):
    """
      Take a stress-controlled step

      Parameters:
        s_np1       next stress
        t_np1       next time
        T_np1       next temperature
    """
    def RJ(e):
      s, h, A = self.model.update_sd(e, np.copy(self.strain_int[-1]),
        T_np1, self.T_int[-1], t_np1, self.t_int[-1], 
        np.copy(self.stress_int[-1]),
        self.stored_int[-1])
      R = s - s_np1
      return R, A

    e_np1 = newton(RJ, np.copy(self.strain_int[-1]), verbose = self.verbose,
        rtol = self.rtol, atol = self.atol, miter = self.miter)

    self.strain_step(np.copy(e_np1), t_np1, T_np1)

  def rate_step(self, sdir, erate, t_np1, T_np1):
    """
      Drive in a given stress direction at a prescribed strain rate, like
      an actual "stress controlled" experiment.

      Parameters:
        sdir        stress direction
        erate       strain rate (in the direction)
        t_np1       next time
        T_np1       next temperature
    """
    sdir /= la.norm(sdir)
    dt = t_np1 - self.t_int[-1]
    def RJ(x):
      a = x[0]
      e_inc = x[1:]
      s, h, A = self.model.update_sd(self.strain_int[-1] + e_inc,
          self.strain_int[-1],
          T_np1, self.T_int[-1], t_np1, self.t_int[-1], self.stress_int[-1],
          self.stored_int[-1])

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
    x0[0] = 1.0
    x0[1:] = sdir / 1000.0
    x = newton(RJ, x0, verbose = self.verbose,
        rtol = self.rtol, atol = self.atol, miter = self.miter)
    e_np1 = x[1:]

    self.strain_step(e_np1, t_np1, T_np1)


class MaximumIterations(RuntimeError):
  def __init__(self):
    super(MaximumIterations, self).__init__("Exceeded the maximum allowed iterations!")

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

  while (nR > rtol * nR0) and (nR > atol * nR):
    x -= la.solve(J, R)
    R, J = RJ(x)
    nR = la.norm(R)
    i += 1
    if verbose:
      print("%i\t%e\t%e\t" % (i, nR, nR / nR0))

    if i > miter:
      print("")
      raise MaximumIterations()

  print("")

  return x
