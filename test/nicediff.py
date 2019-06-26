from __future__ import division

import numpy as np

from common import *

from neml.math import tensors
from neml import history

def diff_scalar_symmetric(fn, S0):
  dfn = lambda s: fn(tensors.Symmetric(s))
  
  return tensors.Symmetric(differentiate(dfn, usym(S0.data))[0])

def diff_history_scalar(fn, h0):
  vec = np.copy(np.array(h0))
  
  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return fn(H)

  nd = differentiate(dfn, vec)

  H = h0.deepcopy()
  H.copy_data(nd)

  return H

def diff_history_symmetric(fn, s0):
  dfn = lambda s: np.array(fn(tensors.Symmetric(usym(s))))
  
  res = differentiate(dfn, s0.data)

  return res

def diff_history_history(fn, h0):
  vec = np.copy(np.array(h0))

  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return np.array(fn(H))

  nd = differentiate(dfn, vec)

  return nd
