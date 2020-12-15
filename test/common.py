from __future__ import division

import numpy as np
import itertools
import os.path

mandel = ((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))
mandel_mults = (1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2))

skew_inds = ((1,2),(0,2),(0,1))
skew_mults = (-1.0,1.0,-1.0)

def localize(fname):
  return os.path.join(os.path.dirname(__file__), fname)

def piece_together_fourth(Dp, Wp):
  """
    Take the skew and symmetric parts of a algorithmic tangent and piece them back together
  """
  sym_id = 0.5*(np.einsum('ik,jl', np.eye(3), np.eye(3)) + np.einsum('jk,il', np.eye(3), np.eye(3)))
  skew_id = 0.5*(np.einsum('ik,jl', np.eye(3), np.eye(3)) - np.einsum('jk,il', np.eye(3), np.eye(3)))

  return np.einsum('ijkl,klmn',sym_tensor_part(Dp), sym_id) + np.einsum('ijkl,klmn', 
      skew_tensor_part(Wp), skew_id) 

def sym_tensor_part(C):
  """
     Take a Mandel stiffness in my notation and convert it back to a full tensor
  """
  Ct = np.zeros((3,3,3,3))
  for a in range(6):
    for b in range(6):
      ind_a = itertools.permutations(mandel[a], r=2)
      ind_b = itertools.permutations(mandel[b], r=2)
      ma = mandel_mults[a]
      mb = mandel_mults[b]
      indexes = tuple(ai+bi for ai, bi in itertools.product(ind_a, ind_b))
      for ind in indexes:
        Ct[ind] = C[a,b] / ma*mb

  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          if l < k:
            Ct[i,j,k,l] = 0.0

  return Ct

def skew_tensor_part(C):
  """
    Take a skew stiffness in my notation and convert it back to a full tensor
  """
  Ct = np.zeros((3,3,3,3))
  for a in range(6):
    for b in range(3):
      inds_a = mandel[a]
      inds_b = skew_inds[b]
      mult_a = mandel_mults[a]
      mult_b = skew_mults[b]
      for ord_a in ((0,1),(1,0)):
        for ord_b, f in zip(((0,1),(1,0)), (1,-1)):
          ind = tuple([inds_a[aa] for aa in ord_a] + [inds_b[bb] for bb in ord_b])
          Ct[ind] = C[a,b] *  mult_a*mult_b * f

  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          if i != j:
            Ct[i,j,k,l] /= 2.0
          if l < k: 
            Ct[i,j,k,l] = 0.0

  return Ct

def unroll_fourth(T):
  """
    Unroll a fourth order tensor into a 9x9
  """
  M = np.zeros((9,9))

  for i in range(3):
    for j in range(3):
      a = i*3+j
      for k in range(3):
        for l in range(3):
          b = k*3+l
          M[a,b] = T[i,j,k,l]

  return M

def reroll_fourth(M):
  """
    Undo unroll_fourth
  """
  T = np.zeros((3,3,3,3))
  for a in range(9):
    i = a // 3
    j = a % 3
    for b in range(9):
      k = b // 3
      l = b % 3
      T[i,j,k,l] = M[a,b]

  return T

def ms2ts(C):
  """
    Convert a Mandel notation stiffness matrix to a full stiffness tensor.
  """
  Ct = np.zeros((3,3,3,3))
  for a in range(6):
    for b in range(6):
      ind_a = itertools.permutations(mandel[a], r=2)
      ind_b = itertools.permutations(mandel[b], r=2)
      ma = mandel_mults[a]
      mb = mandel_mults[b]
      indexes = tuple(ai+bi for ai, bi in itertools.product(ind_a, ind_b))
      for ind in indexes:
        Ct[ind] = C[a,b] / (ma*mb)

  return Ct

def ts2ms(C):
  """
    Convert a stiffness tensor into a Mandel notation stiffness matrix
  """
  Cv = np.zeros((6,6))
  for i in range(6):
    for j in range(6):
      ma = mandel_mults[i]
      mb = mandel_mults[j]
      Cv[i,j] = C[mandel[i]+mandel[j]] * ma * mb

  return Cv

def ws2ts(C):
  """
    Convert a skew notation stiffness matrix to a full stiffness tensor.
  """
  Ct = np.zeros((3,3,3,3))
  for a in range(6):
    for b in range(3):
      inds_a = mandel[a]
      inds_b = skew_inds[b]
      mult_a = mandel_mults[a]
      mult_b = skew_mults[b]
      for ord_a in ((0,1),(1,0)):
        for ord_b, f in zip(((0,1),(1,0)), (1,-1)):
          ind = tuple([inds_a[aa] for aa in ord_a] + [inds_b[bb] for bb in ord_b])
          Ct[ind] = C[a,b] / (mult_a*mult_b) * f

  return Ct

def ts2ws(C):
  """
    Convert a stiffness tensor into a skew notation stiffness matrix
  """
  Cv = np.zeros((6,3))
  for i in range(6):
    for j in range(3):
      ma = mandel_mults[i]
      mb = skew_mults[j]
      Cv[i,j] = C[mandel[i]+skew_inds[j]] * ma * mb

  return Cv

def wws2ts(C):
  """
    Convert a d_skew_d_sym matrix to a tensor
  """
  Ct = np.zeros((3,3,3,3))
  for a in range(3):
    for b in range(6):
      inds_a = skew_inds[a]
      inds_b = mandel[b]
      mult_a = skew_mults[a]
      mult_b = mandel_mults[b]
      for ord_a, f in zip(((0,1),(1,0)),(1,-1)):
        for ord_b in ((0,1),(1,0)):
          ind = tuple([inds_a[aa] for aa in ord_a] + [inds_b[bb] for bb in ord_b])
          Ct[ind] = C[a,b] / (mult_a * mult_b) * f

  return Ct

def ts2wws(C):
  """
    Convert a stiffness tensor into a d_skew_d_sym matrix
  """
  Cv = np.zeros((3,6))
  for i in range(3):
    for j in range(6):
      ma = skew_mults[i]
      mb = mandel_mults[j]
      Cv[i,j] = C[skew_inds[i] + mandel[j]] * ma * mb

  return Cv

def ts2sww(C):
  """
    Convert a stiffness tensor into a d_sym_d_skew matrix
  """
  Cv = np.zeros((6,3))
  for i in range(6):
    for j in range(3):
      ma = mandel_mults[i]
      mb = skew_mults[j]
      Cv[i,j] = C[mandel[i] + skew_inds[j]] * ma * mb

  return Cv

def sym(A):
  """
    Take a symmetric matrix to the Mandel convention vector.
  """
  return np.array([A[0,0], A[1,1], A[2,2], np.sqrt(2)*A[1,2], 
    np.sqrt(2)*A[0,2], np.sqrt(2)*A[0,1]])

def usym(v):
  """
    Take a Mandel symmetric vector to the full matrix.
  """
  return np.array([
    [v[0], v[5]/np.sqrt(2), v[4]/np.sqrt(2)],
    [v[5]/np.sqrt(2), v[1], v[3]/np.sqrt(2)],
    [v[4]/np.sqrt(2), v[3]/np.sqrt(2), v[2]]
    ])

def skew(A):
  """
    Take a skew matrix to my vector convention.
  """
  return np.array([-A[1,2], A[0,2], -A[0,1]])

def uskew(v):
  """
    Take a skew vector in my convention to the matrix
  """
  return np.array([
    [0, -v[2], v[1]],
    [v[2], 0, -v[0]],
    [-v[1], v[0], 0]
    ])

def differentiate(fn, x0, eps = 1.5e-6):
  if np.isscalar(x0):
    sx = (1,)
  else:
    sx = x0.shape

  f0 = fn(x0)

  if np.isscalar(f0):
    sf = (1,)
  else:
    sf = f0.shape

  tshape = sf + sx
  rshape = sx + sf

  D = np.zeros(rshape)

  for index in np.ndindex(sx):
    if np.isscalar(x0):
      diff = np.abs(x0 * eps)
    else:
      diff = np.abs(x0[index] * eps)
    if diff < eps:
      diff = eps
    
    if np.isscalar(x0):
      xp = diff
    else:
      xp = np.zeros(sx)
      xp[index] = diff

    fp = fn(xp+x0)

    D[index] = (fp - f0) / diff
  
  Df = np.zeros(tshape)
  # Reverse
  for ind1 in np.ndindex(sx):
    for ind2 in np.ndindex(sf):
      Df[ind2+ind1] = D[ind1+ind2]
  
  if np.isscalar(x0):
    return Df[0][0]
  else:
    return Df

def make_dev(s):
  return s - np.array([1,1,1,0,0,0]) * np.sum(s[:3]) / 3.0
