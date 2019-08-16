import numpy as np

from neml.math import rotations
from neml.cp import harmonics

import numpy.linalg as la
import scipy.optimize as opt
import scipy.special as spfns

def hiter(order):
  return ((n,i,j) for n in range(order+1) 
      for i in range(-n,n+1) for j in range(-n,n+1))

class ODF(object):
  """
    Parent class of all orientation distribution function implementations
  """
  def __init__(self):
    pass

  @property
  def pdf2mrd(self):
    return 8 * np.pi**2.0

class HarmonicsODF(ODF):
  """
    Orientation distribution function based on the generalized spherical
    harmonics
  """
  def __init__(self, order, *args, **kwargs):
    """
      Parameters:
        order       order of the harmonic expansion to use
    """
    super().__init__(*args, **kwargs)
    self.order = order
    self.p = np.zeros((harmonics.nharm_total(self.order),), dtype = complex)

  @property
  def ncoefs(self):
    return len(self.p)

  def value(self, pt):
    return np.real(sum(harmonics.harmonic_SO3(n,i,j,pt) * self.p[it] for 
        it,(n,i,j) in enumerate(hiter(self.order))))

  def project(self, pts, wts = None, method = 'nnls'):
    """
      Project discrete points onto the harmonics

      Parameters:
        pts         points, as quaternions

      Optional:
        wts         point weights, defaults to all ones
        method      which method to use for projection
                    Options:
                      series        the standard Fourier method
                      ls            least squares
                      nnls          non-negative least squares
    """
    if wts is None:
      wts = np.ones((len(pts),))

    if len(pts) != len(wts):
      raise ValueError("Length of pts and wts should be the same")
    
    if method == 'series':
      self.project_series(pts, wts)
    elif method == 'ls':
      self.project_ls(pts, wts)
    elif method == 'nnls':
      self.project_nnls(pts, wts)
    else:
      raise ValueError("Unknown projection method %s!" % method)

  def project_series(self, pts, wts):
    """
      Project using the standard Dirac delta/Fourier method
    """
    # Project
    for wt,pt in zip(wts,pts):
      for it,(n,i,j) in enumerate(hiter(self.order)):
        self.p[it] += wt * np.conj(harmonics.harmonic_SO3(n,i,j,pt))
    
    # Normalize
    self.p /= np.sum(wts)

  def project_ls(self, pts, wts):
    """
      Project using least squares.
    """
    A = np.array([[harmonics.harmonic_SO3(n,i,j,pt) for n,i,j in 
      hiter(self.order)] for pt in pts])
    self.p = la.lstsq(A, wts, rcond = None)[0]

    # Normalize
    qpts, qwts = harmonics.quadrature_SO3(self.order)
    v = sum(self.value(pt) * wt for wt,pt in zip(qwts,qpts))
    self.p /= v

  def project_nnls(self, pts, wts):
    """
      Project using nonnegative least squares
    """
    Ap = np.array([[harmonics.harmonic_SO3(n,i,j,pt) for n,i,j in 
      hiter(self.order)] for pt in pts])
    A = np.vstack((np.real(Ap),np.imag(Ap)))
    b = np.hstack((wts, np.zeros(wts.shape)))

    self.p = opt.nnls(A,b)[0]

    # Normalize
    qpts, qwts = harmonics.quadrature_SO3(self.order)
    v = sum(self.value(pt) * wt for wt,pt in zip(qwts,qpts))
    self.p /= v

  def pole_figure_harmonics(self, poles):
    """
      Return the coefficients for the spherical harmonics for the given
      poles up to the degree of the 

      Parameters
        pole        poles in a nx3 array as real vectors in the right 
                    coordinates
    """
    pass
