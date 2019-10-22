from __future__ import division

import numpy as np

from common import *

from neml.math import tensors
from neml import history

def diff_skew_symmetric(fn, s0):
  dfn = lambda s: fn(tensors.Symmetric(usym(s))).data

  return tensors.SkewSymR4(differentiate(dfn, s0.data))

def diff_skew_history(fn, h0):
  vec = np.copy(np.array(h0))

  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return fn(H).data

  return differentiate(dfn, vec)

def diff_symmetric_symmetric(fn, s0):
  dfn = lambda s: fn(tensors.Symmetric(usym(s))).data

  return tensors.SymSymR4(differentiate(dfn, s0.data))

def diff_symmetric_skew(fn, w0):
  dfn = lambda w: fn(tensors.Skew(uskew(w))).data

  return tensors.SymSkewR4(differentiate(dfn, w0.data))

def diff_symmetric_history(fn, h0):
  vec = np.copy(np.array(h0))

  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return fn(H).data

  nd = differentiate(dfn, vec)

  return nd

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

def diff_history_skew(fn, w0):
  dfn = lambda w: np.array(fn(tensors.Skew(uskew(w))))

  return differentiate(dfn, w0.data)

def diff_history_history(fn, h0):
  vec = np.copy(np.array(h0))

  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return np.array(fn(H))

  nd = differentiate(dfn, vec)

  return nd
