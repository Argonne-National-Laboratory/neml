#!/usr/bin/env python

from __future__ import division

import numpy as np
import copy
import itertools

from sympy import *

mandel = ((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))
mandel_mults = (1,1,1,sqrt(2),sqrt(2),sqrt(2))

skew_inds = ((1,2),(0,2),(0,1))
skew_mults = (-1,1,-1)

def object_einsum(string, *arrays):
  """Simplified object einsum, not as much error checking
    
  does not support "..." or list input and will see "...", etc. as three times
  an axes identifier, tries normal einsum first!
    
  NOTE: This is untested, and not fast, but object type is
  never really fast anyway...
  """
  try:
    return np.einsum(string, *arrays)
  except TypeError:
    pass
  
  s = string.split('->')
  in_op = s[0].split(',')
  out_op = None if len(s) == 1 else s[1].replace(' ', '')

  in_op = [axes.replace(' ', '') for axes in in_op]
  all_axes = set()

  for axes in in_op:
    all_axes.update(axes)

  if out_op is None:
    out_op = sorted(all_axes)
  else:
    all_axes.update(out_op)
  
  perm_dict = {_[1]: _[0] for _ in enumerate(all_axes)}
  
  dims = len(perm_dict)
  op_axes = []
  for axes in (in_op + list((out_op,))):
    op = [-1] * dims
    for i, ax in enumerate(axes):
      op[perm_dict[ax]] = i
    op_axes.append(op)
  
  op_flags = [('readonly',)] * len(in_op) + [('readwrite', 'allocate')]
  dtypes = [np.object_] * (len(in_op) + 1) # cast all to object

  nditer = np.nditer(arrays + (None,), op_axes=op_axes, flags=['buffered', 'delay_bufalloc', 'reduce_ok', 'grow_inner', 'refs_ok'], op_dtypes=dtypes, op_flags=op_flags)

  nditer.operands[-1][...] = 0
  nditer.reset()
  
  for vals in nditer:
    out = vals[-1]
    prod = copy.deepcopy(vals[0])
    for value in vals[1:-1]:
      prod *= value
    out += prod
    
  return nditer.operands[-1]

def zero_tensor(shape):
  l = np.prod(shape)
  return MutableDenseNDimArray([0]*l, shape)

def piece_together_fourth(Dp, Wp):
  """
    Take the skew and symmetric parts of a algorithmic tangent and piece them back together
  """
  sym_id = (object_einsum('ik,jl', eye(3), eye(3)) + object_einsum('jk,il', eye(3), eye(3)))/2
  skew_id = (object_einsum('ik,jl', eye(3), eye(3)) - object_einsum('jk,il', eye(3), eye(3)))/2

  D = sym_tensor_part(Dp)
  W = skew_tensor_part(Wp)
  res = zero_tensor((3,3,3,3))

  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          for m in range(3):
            for n in range(3):
              res[i,j,m,n] += D[i,j,k,l] * sym_id[k,l,m,n] + W[i,j,k,l] * skew_id[k,l,m,n]

  return res

def sym_tensor_part(C):
  """
     Take a Mandel stiffness in my notation and convert it back to a full tensor
  """
  Ct = zero_tensor((3,3,3,3))
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
            Ct[i,j,k,l] = 0

  return Ct

def skew_tensor_part(C):
  """
    Take a skew stiffness in my notation and convert it back to a full tensor
  """
  Ct = zero_tensor((3,3,3,3))
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
            Ct[i,j,k,l] /= 2
          if l < k: 
            Ct[i,j,k,l] = 0

  return Ct

def unroll_fourth(T):
  M = zeros(9,9)

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
  T = zero_tensor((3,3,3,3))
  for a in range(9):
    i = a // 3
    j = a % 3
    for b in range(9):
      k = b // 3
      l = b % 3
      T[i,j,k,l] = M[a,b]

  return T

def unroll_second(T):
  return Matrix([T[i,j] for i in range(3) for j in range(3)])

def ms2ts(C):
  """
    Convert a Mandel notation stiffness matrix to a full stiffness tensor.
  """
  Ct = zero_tensor((3,3,3,3))
  for a in range(6):
    for b in range(6):
      ind_a = itertools.permutations(mandel[a], r=2)
      ind_b = itertools.permutations(mandel[b], r=2)
      ma = mandel_mults[a]
      mb = mandel_mults[b]
      indexes = tuple(ai+bi for ai, bi in itertools.product(ind_a, ind_b))
      for i,j,k,l in indexes:
        Ct[i,j,k,l] = C[a,b] / (ma*mb)

  return Ct

def ts2ms(C):
  """
    Convert a stiffness tensor into a Mandel notation stiffness matrix
  """
  Cv = zeros(6,6)
  for i in range(6):
    for j in range(6):
      ma = mandel_mults[i]
      mb = mandel_mults[j]
      Cv[i,j] = C[mandel[i]+mandel[j]] * ma * mb

  return Cv

def ts2ms2(C):
  """
    Convert a stiffness tensor into a Mandel notation stiffness matrix
  """
  Cv = zeros(6,6)
  for i in range(6):
    for j in range(6):
      ma = mandel_mults[i]
      mb = mandel_mults[j]
      Cv[i,j] = (C[mandel[i]+mandel[j]] + C[mandel[i]+mandel[j][::-1]] + C[mandel[i][::-1]+mandel[j]] + C[mandel[i][::-1]+mandel[j][::-1]]) * ma * mb / 4

  return Cv

def ws2ts(C):
  """
    Convert a skew notation stiffness matrix to a full stiffness tensor.
  """
  Ct = zero_tensor((3,3,3,3))
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
  Cv = zeros(6,3)
  for i in range(6):
    for j in range(3):
      ma = mandel_mults[i]
      mb = skew_mults[j]
      Cv[i,j] = C[mandel[i]+skew_inds[j]] * ma * mb

  return Cv

def wws2ts(C):
  """
    Convert a skew notation stiffness matrix to a full stiffness tensor.
  """
  Ct = zero_tensor((3,3,3,3))
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
    Convert a stiffness tensor into a skew notation stiffness matrix
  """
  Cv = zeros(3,6)
  for i in range(3):
    for j in range(6):
      ma = skew_mults[i]
      mb = mandel_mults[j]
      Cv[i,j] = C[skew_inds[i] + mandel[j]] * ma * mb

  return Cv

def ts2wsw(C):
  """
    Convert a stiffness tensor into a skew thing
  """
  Cv = zeros(6,3)
  for i in range(6):
    for j in range(3):
      ma = mandel_mults[i]
      mb = skew_mults[j]
      Cv[i,j] = C[mandel[i] + skew_inds[j]] * ma * mb

  return Cv

def trace(X):
  return (X[0,0] + X[1,1] + X[2,2])

def sym(A):
  return [A[0,0], A[1,1], A[2,2], sqrt(2) * A[1,2], sqrt(2) * A[0,2], sqrt(2) * A[0,1]]

if __name__ == "__main__":
  I = Matrix([[1,0,0],[0,1,0],[0,0,1]])

  n1, n2, n3 = symbols("n[0] n[1] n[2]")
  n = Matrix([n1,n2,n3])

  def sack(expr):
    return expr.replace(
        lambda x: x.is_Pow and x.exp > 0 and not x.base.is_integer,
        lambda x: Symbol('*'.join([x.base.name]*x.exp)))

  normal = zero_tensor((3,3,3,3))

  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          normal[i,j,k,l] = n[k] * n[l] * n[i] * n[j] 

  normal = simplify(normal)

  normal_sym = zero_tensor((3,3,3,3))
  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          normal_sym[i,j,k,l] = (normal[i,j,k,l] + normal[j,i,k,l])/2
  normal_sym = simplify(ts2ms2(normal_sym))
  
  print("First operator")
  for i in range(6):
    for j in range(6):
      print(("\tSS[%i,%i] = " + str(sack(normal_sym[i,j])) + ";") % (i,j))     


  shear = zero_tensor((3,3,3,3))
  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          shear[i,j,k,l] = (I[i,k] - n[i]*n[k]) * n[j] * n[l]
  shear = simplify(shear)

  shear_sym = zero_tensor((3,3,3,3))
  for i in range(3):
    for j in range(3):
      for k in range(3):
        for l in range(3):
          shear_sym[i,j,k,l] = (shear[i,j,k,l] + shear[j,i,k,l])

  shear_sym = simplify(ts2ms2(shear_sym))
  
  print("Second operator")
  for i in range(6):
    for j in range(6):
      print(("\tSS[%i,%i] = " + str(sack(shear_sym[i,j])) + ";") % (i,j))     

