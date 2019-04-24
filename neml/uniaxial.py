import numpy as np
import numpy.linalg as la
import scipy.optimize as opt
from neml.nlsolvers import newton

class UniaxialModel(object):
  """
    Takes a NEML model as input and performs the necessary calculations
    to give it a uniaxial response.
  """
  def __init__(self, nemlmodel, verbose = False, rtol = 1.0e-6, atol = 1.0e-10,
      miter = 20):
    """
      Parameters:
        nemlmodel:       the underlying 3D NEML model

      Keyword Args:
        verbose:         verbosity flag
        rtol:            newton relative tolerance
        atol:            newton absolute tolerance
        miter:           newton max iterations
    """
    self.model = nemlmodel
    self.verbose = verbose
    self.rtol = rtol
    self.atol = atol
    self.miter = miter

  @property
  def nstore():
    """
      Number of variables required to support the material update
    """
    return self.model.nstore + 5

  def init_store(self):
    """
      Initialize the variables required to support the material update
    """
    return np.concatenate((np.zeros((5,)), self.model.init_store()))

  @property
  def nhist():
    """
      True number of model history variables
    """
    return self.model.nhist

  def init_hist(self):
    """
      Initialize the model history
    """
    return self.model.init_hist()

  def alpha(self, T):
    """
      Returns the temperature-dependent thermal expansion coefficient

      Parameters:
        T       temperature
    """
    return self.model.alpha(T)

  def elastic_strains(self, s, T, h):
    """
      Returns the elastic strain for a given stress and temperature
    """
    s_vec = self.ss2sv(s)
    return self.model.elastic_strains(s_vec, T, h)[0]

  def update(self, e_np1, e_n, T_np1, T_n, t_np1, t_n, s_n, h_n, u_n, p_n):
    """
      Material update function
      
      Parameters:
        e_np1:       next strain
        e_n:         previous strain
        T_np1:       next temperature
        T_n:         previous temperature
        t_np1:       next time
        t_n:         previous time
        s_n:         previous stress
        h_n:         previous history
        u_n:         previous energy
        p_n:         previous dissipation
    """
    e_n_vec = self.es2ev(e_n, h_n)
    s_n_vec = self.ss2sv(s_n)
    h_n_vec = h_n[5:]
    
    def RJ(x):
      e_np1_vec = self.es2ev(e_np1, x)
      s_np1, h_np1, A_np1, u_np1, p_np1 = self.model.update_sd(
          e_np1_vec, e_n_vec, T_np1, T_n, t_np1, t_n, s_n_vec, h_n_vec,
          u_n, p_n)
      
      R = s_np1[1:]
      J = A_np1[1:,:][:,1:]

      return R, J

    x0 = h_n[:5]
    ef = newton(RJ, x0, verbose = self.verbose, rtol = self.rtol,
        atol = self.atol, miter = self.miter)
    
    e_np1_vec = self.es2ev(e_np1, ef)
    s_np1_vec, h_np1_vec, A_np1_vec, u_np1, p_np1 = self.model.update_sd(
        e_np1_vec, e_n_vec, T_np1, T_n, t_np1, t_n, s_n_vec, h_n_vec,
        u_n, p_n)
    h_np1 = np.concatenate((ef, h_np1_vec))
    s_np1 = s_np1_vec[0]

    # Faster than inverting
    sv = la.solve(A_np1_vec[1:,:][:,1:], A_np1_vec[1:,0])
    A_np1 = A_np1_vec[0,0] - np.dot(A_np1_vec[0,1:], sv)

    return s_np1, h_np1, A_np1, u_np1, p_np1

  def ss2sv(self, s):
    """
      Helper to convert scalar stress into stress vector

      Parameters:
        s:       scalar stress
    """
    vec = np.zeros((6,))
    vec[0] = s
    return vec

  def es2ev(self, e, h):
    """
      Helper to convert scalar strains into strain vectors

      Parameters:
        e:       scalar strain
        h:       history (stores other 5 strain components)
    """
    vec = np.zeros((6,))
    vec[0] = e
    vec[1:] = h[:5]
    return vec
