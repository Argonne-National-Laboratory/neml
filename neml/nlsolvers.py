import numpy as np
import numpy.linalg as la

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

def newton_precalc(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-10, 
    miter = 50):
  """
    Manually-coded newton-raphson with precalculated jacobian

    This version doesn't do linesearch because ironically that would 
    require re-inverting everything.

    Parameters:
      RJ            function return the residual + a function that gives
                    the action of the inverse jacobian on a vector
      x0            initial guess

    Optional:
      verbose       verbose output
      rtol          relative tolerance
      atol          absolute tolerance
      miter         maximum iterations
  """
  R, J = RJ(x0)
  nR = la.norm(R)
  nR0 = nR
  x = np.copy(x0)
  
  i = 0

  if verbose:
    print("Iter.\tnR\t\tnR/nR0")
    print("%i\t%e\t%e\t" % (i, nR, 1.0))

  while (nR > rtol * nR0) and (nR > atol):
    a = J(R)
    x -= a
    R, J = RJ(x)
    nR = la.norm(R)
    i += 1
    if verbose:
      print("%i\t%e\t%e" % (i, nR, nR / nR0))
    if i > miter:
      if verbose:
        print("")
      raise MaximumIterations()
  
  if verbose:
    print("")

  return x

def newton(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-10, miter = 50,
    linesearch = 'none', bt_tau = 0.5, bt_c = 1.0e-4):
  """
    Manually-code newton-raphson so that I can output convergence info, if
    requested.

    Parameters:
      RJ            function return the residual + jacobian
      x0            initial guess

    Optional:
      verbose       verbose output
      rtol          relative tolerance
      atol          absolute tolerance
      miter         maximum iterations
      linesearch    available options: "none" and "backtracking"
      bt_tau        tau factor for backtracking line search
      bt_c          c factor for backtracking line search
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
      rtol      relative tolerance
      atol      absolute tolerance
      miter     maximum number of iterations
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

def scalar_newton(RJ, x0, verbose = False, rtol = 1.0e-6, atol = 1.0e-10,
    miter = 20, quasi = None):
  """
    Manually-code newton-raphson for scalar problems

    Parameters:
      RJ        function return the residual + jacobian
      x0        initial guess

    Optional:
      verbose   verbose output
      rtol      relative tolerance
      atol      absolute tolerance
      miter     maximum solver iterations
      quasi     if true use a quasi-newton method
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

