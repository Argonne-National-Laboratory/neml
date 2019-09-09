#!/usr/bin/env python3

import numpy as np
import numpy.polynomial as npp
from math import factorial

from sympy import N
from sympy.physics.quantum.spin import Rotation

def rs(i,l):
  return i + l

def wigner_d_1(l, x):
  dvals = []
  for ll in range(l+1):
    dvals.append(np.zeros((2*ll+1,2*ll+1)))
    
    for k in range(-ll,ll+1):
      dvals[ll][rs(0,ll),rs(k,ll)] = 

  return dvals

def P_l(l,x):
  pvals = [np.array([1.0])]
  for ll in range(1,l+1):
    pvals.append(np.zeros((ll+1,)))
    pvals[ll][ll] = np.sqrt(float(2*ll-1)/(2*ll))*np.sin(x) * pvals[ll-1][ll-1]
    pvals[ll][ll-1] = np.sqrt(float(2*ll-1))*np.cos(x) * pvals[ll-1][ll-1]
    for k in range(0,ll-1):
      pvals[ll][k] = ((2.0*ll-1)*np.cos(x) * pvals[ll-1][k] - 
          np.sqrt(float(ll-k-1) *(ll+k-1)) * pvals[ll-2][k]) / np.sqrt(float(ll-k)*(ll+k))

  return pvals

def P_1_recurse(l, k, x):
  if l == k == 0:
    return 1.0
  elif l == k:
    return np.sqrt((2.0*l-1.0)/(2.0*l)) * np.sin(x) * P_1_recurse(l-1,l-1, x)
  elif k == (l-1):
    return np.sqrt(2.0*l - 1.0) * np.cos(x) * P_1_recurse(l-1, k, x)
  else:
    return ((2.0*l-1.0) * np.cos(x) * P_1_recurse(l-1,k, x) 
        - np.sqrt((l-k-1.0)*(l+k-1.0)) * P_1_recurse(l-2,k,x)) / np.sqrt(float((l-k)*(l+k)))

def reference(l, m, k, x):
  return complex(N(Rotation.d(l, m, k, x).doit()))

if __name__ == "__main__":
  #print(reference(0,0,0,np.pi/4))
  #print(P_1_recurse(2,2, np.pi/8.0))
  #print(P_l(2, np.pi/8.0))

  print(wigner_d_1(2, np.pi/8.0))
