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

def trace(X):
  return (X[0,0] + X[1,1] + X[2,2])

def sym(A):
  return [A[0,0], A[1,1], A[2,2], sqrt(2) * A[1,2], sqrt(2) * A[0,2], sqrt(2) * A[0,1]]

def update():
  W0, W1, W2 = symbols("W[0] W[1] W[2]")
  D0, D1, D2, D3, D4, D5 = symbols("D[0] D[1] D[2] D[3] D[4] D[5]")
  S0, S1, S2, S3, S4, S5 = symbols("Sn[0] Sn[1] Sn[2] Sn[3] Sn[4] Sn[5]")
  O0, O1, O2, O3, O4, O5 = symbols("So[0] So[1] So[2] So[3] So[4] So[5]")

  G00, G01, G02, G10, G11, G12, G20, G21, G22 = symbols("G00 G01 G02 G10 G11 G12 G20 G21 G22")

  W = Matrix([[0,-W2,W1],[W2,0,-W0],[-W1,W0,0]])
  D = Matrix([[D0, D5/sqrt(2), D4/sqrt(2)],[D5/sqrt(2),D1,D3/sqrt(2)],[D4/sqrt(2),D3/sqrt(2),D2]])
  L = D + W

  S = Matrix([[S0, S5/sqrt(2), S4/sqrt(2)],[S5/sqrt(2),S1,S3/sqrt(2)],[S4/sqrt(2),S3/sqrt(2),S2]])
  O = Matrix([[O0, O5/sqrt(2), O4/sqrt(2)],[O5/sqrt(2),O1,O3/sqrt(2)],[O4/sqrt(2),O3/sqrt(2),O2]])

  G = Matrix([[G00,G01,G02],[G10,G11,G12],[G20,G21,G22]])
  
  V1 = simplify(O + S*(L.T) + L*S - (trace(L) * S))

  V11 = simplify(sym(V1))
  
  print("RHS")
  for i in range(6):
    print(("\tSt[%i] = " + str(V11[i]) + ";") % i)
  print("")

  J = object_einsum('im,jn', eye(3), eye(3)) * (1 + trace(L)) - object_einsum('im,jn', eye(3), L) - object_einsum('im,jn', L, eye(3))
  JM = simplify(unroll_fourth(J))
  
  print("Matrix")
  for i in range(9):
    for j in range(9):
      print(("\tM[%i] = " + str(JM[i,j]) + ";") % (i*9+j))
  print("")

def tangent():
  sym_mat = Matrix([[symbols("M[%i]" % (i*6+j)) for j in range(6)] for i in range(6)])
  sym_ten = ms2ts(sym_mat)
  sym_full = unroll_fourth(sym_ten)

  skew_mat = Matrix([[symbols("M[%i]" % (i*3+j)) for j in range(3)] for i in range(6)])
  skew_ten = ws2ts(skew_mat)
  skew_full = unroll_fourth(skew_ten)

  A_mat = Matrix([[symbols("A[%i]" % (i*9+j)) for j in range(9)] for i in range(9)])
  A_ten = reroll_fourth(A_mat)
  A_mandel = ts2ms(A_ten)
  A_skew = ts2ws(A_ten)

  print("Mandel->9x9")
  for i in range(9):
    for j in range(9):
      print(("\tA[%i] = " + str(sym_full[i,j]) + ";") % (i*9+j))
  print("")
  print("9x9->Mandel")
  for i in range(6):
    for j in range(6):
      print(("\tM[%i] = " + str(A_mandel[i,j]) + ";") % (i*6+j))
  print("")
  
  print("Skew->9x9")
  for i in range(9):
    for j in range(9):
      print(("\tA[%i] = " + str(skew_full[i,j]) + ";") % (i*9+j))
  print("")
  
  print("9x9->Skew")
  for i in range(6):
    for j in range(3):
      print(("\tM[%i] = " + str(A_skew[i,j]) + ";") % (i*3+j))
  print("")

  S0, S1, S2, S3, S4, S5 = symbols("S[0] S[1] S[2] S[3] S[4] S[5]")
  S = Matrix([[S0, S5/sqrt(2), S4/sqrt(2)],[S5/sqrt(2),S1,S3/sqrt(2)],[S4/sqrt(2),S3/sqrt(2),S2]])
  tensor = simplify(object_einsum('in,jm', S, eye(3)) + object_einsum('im,nj', eye(3), S) - object_einsum('ij,mn', S, eye(3)))
  matrix = simplify(unroll_fourth(tensor))
  print("Matrix")
  for i in range(9):
    for j in range(9):
      print(("\tM[%i] = " + str(matrix[i,j]) + ";") % (i*9+j))
  print("")

def parts_to_whole():
  sym_mat = Matrix([[symbols("D[%i]" % (i*6+j)) for j in range(6)] for i in range(6)])
  skew_mat = Matrix([[symbols("W[%i]" % (i*3+j)) for j in range(3)] for i in range(6)])
  
  tensor = piece_together_fourth(sym_mat, skew_mat)
  matrix = simplify(unroll_fourth(tensor))
  print("Matrix")
  for i in range(9):
    for j in range(9):
      print(("\tM[%i] = " + str(matrix[i,j]) + ";") % (i*9+j))
  print("")

if __name__ == "__main__":
  init_printing(use_unicode=True)

  update()
  tangent()
  parts_to_whole()


