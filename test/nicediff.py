from __future__ import division

import numpy as np

from common import *

from neml.math import tensors
from neml import history

def diff_symsym_sym(fn, s0):
    dfn = lambda s: fn(tensors.Symmetric(usym(s))).data

    return tensors.SymSymSymR6(differentiate_new(dfn, s0.data).reshape(6,6,6))

def diff_skew_symmetric(fn, s0):
  dfn = lambda s: fn(tensors.Symmetric(usym(s))).data

  return tensors.SkewSymR4(differentiate_new(dfn, s0.data))

def diff_skew_history(fn, h0):
  vec = np.copy(np.array(h0))

  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return fn(H).data

  return differentiate_new(dfn, vec)

def diff_symmetric_symmetric(fn, s0):
  dfn = lambda s: fn(tensors.Symmetric(usym(s))).data

  return tensors.SymSymR4(differentiate_new(dfn, s0.data))

def diff_symmetric_skew(fn, w0):
  dfn = lambda w: fn(tensors.Skew(uskew(w))).data

  return tensors.SymSkewR4(differentiate_new(dfn, w0.data))

def diff_symsym_history(fn, h0):
  vec = np.copy(np.array(h0))

  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return fn(H).data

  nd = differentiate_new(dfn, vec)

  return nd

def diff_symmetric_history(fn, h0):
  vec = np.copy(np.array(h0))

  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return fn(H).data

  nd = differentiate_new(dfn, vec)

  return nd

def diff_symmetric_scalar(fn, s0):
  dfn = lambda s: fn(s).data

  res = differentiate_new(dfn, np.array([s0]))

  return tensors.Symmetric(usym(res[:,0]))

def diff_scalar_symmetric(fn, S0):
  dfn = lambda s: fn(tensors.Symmetric(s))
  
  return tensors.Symmetric(differentiate_new(dfn, usym(S0.data))[0])

def diff_history_scalar(fn, h0):
  vec = np.copy(np.array(h0))
  
  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return fn(H)

  nd = differentiate_new(dfn, vec)

  H = h0.deepcopy()
  H.copy_data(nd)

  return H

def diff_history_symmetric(fn, s0):
  dfn = lambda s: np.array(fn(tensors.Symmetric(usym(s))))
  
  res = differentiate_new(dfn, s0.data)

  return res

def diff_history_skew(fn, w0):
  dfn = lambda w: np.array(fn(tensors.Skew(uskew(w))))

  return differentiate_new(dfn, w0.data)

def diff_history_history(fn, h0):
  vec = np.copy(np.array(h0))

  def dfn(x):
    H = h0.deepcopy()
    H.copy_data(x)
    return np.array(fn(H))

  nd = differentiate_new(dfn, vec)

  return nd
