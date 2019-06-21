from __future__ import division

import numpy as np

from common import *

from neml.math import tensors
from neml import history

def diff_scalar_symmetric(fn, S0):
  dfn = lambda s: fn(tensors.Symmetric(s))
  
  return tensors.Symmetric(differentiate(dfn, usym(S0.data))[0])

def diff_scalar_history(fn, h0):
  vec = np.copy(np.array(h0))
  
  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return fn(H)

  nd = differentiate(dfn, vec)

  H = h0.deepcopy()
  H.copy_data(nd)

  return H
