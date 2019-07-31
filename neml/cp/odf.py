import numpy as np

from neml.math import rotations

class ODF(object):
  """
    Parent class of all orientation distribution function implementations
  """
  def __init__(self):
    pass

class HarmonicsODF(object):
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

  def project(self, pts, wts = None):
    """
      Project discrete points onto the harmonics

      Parameters:
        pts         points, as quaternions

      Optional:
        wts         point weights, defaults to all ones
    """
    if wts is None:
      wts = np.ones((len(pts),))

    if len(pts) != len(wts):
      raise ValueError("Length of pts and wts should be the same")


