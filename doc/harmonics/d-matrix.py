#!/usr/bin/env python3

import numpy as np
import numpy.polynomial as npp
from math import factorial

from spherical_functions import Wigner_D_element

def pder(l,m,n):
  """
    The derivative in the Legendre polynomials

    Parameters:
      l,m,n     order
  """
  roots = [1.0] * (l - m) + [-1.0] * (l + m)
  poly = npp.polynomial.polyfromroots(roots)
  pder = npp.polynomial.polyder(poly, m = (l - n))

  return pder

def P_x(x, l, m, n):
  """
    The associated Legendre polynomials, in Bunge's convention

    Parameters:
      x         variable
      l,m,n     order
  """
  return ((-1.0)**(l-m) * (1.0j)**(n-m) / (float(factorial(l-m))*2.0**l)
      ) * np.sqrt(float(factorial(l-m)*factorial(l+n))/float(factorial(l+m)*factorial(l-n))
          ) * ((1.0 - x)**(-(n-m)/2.0) * (1.0 + x)**(-(n+m)/2.0)
              ) * npp.polynomial.polyval(x, pder(l,m,n)) 

def P(t, l, m, n):
  return P_x(np.cos(t), l, m, n)

def P_alt(t, l, m, n):
  return Wigner_D_element(0,t,0,l,m,n) / ((-1.0)**(l-m) * (1.0j)**(n-m))

def W_all(t, lmax, eps = np.finfo(float).eps / 2):
  """
    Compute all the Wigner d matrices up to order l 
  """
  pass

def W_comp_all(t, lmax):
  """
    Actual calculation function for 0 < t < pi/2
  """
  # Get the g values
  g = [np.zeros((l+1,)) for l in range(lmax+1)]
  g[0][0] = 1.0
  for l in range(1,lmax+1):
    g[l][0] = np.sqrt(float(2*l-1)/(2*l)) * g[l-1][0]
    for m in range(1,l+1):
      g[l][m] = np.sqrt(float(l-m+1)/(l+m)) * g[l][m-1]
  
  zos1 = lambda i, l: i + l
  zos2 = lambda i, l: i + l
  vals = [np.zeros((2*i+1,2*i+1), dtype = complex) for i in range(lmax+1)]
  
  # Precalc some cosines and sines
  c = np.cos(t)
  s = np.sin(t)

  # n = l with (28)
  for l in range(0,lmax+1):
    for m in range(0,l+1):
      vals[l][zos1(m,l),zos2(l,l)] = (-1)**(l+m) * g[l][m] * (1.0+c)**m * s**(l-m)

  # m = l with (26)
  for l in range(0,lmax+1):
    for n in range(l,-l,-1):
      vals[l][zos1(l,l),zos2(n-1,l)] = (l+n)/np.sqrt(float(l*(l+1)-n*(n-1))) * s / (1+c) * vals[l][zos1(l,l),zos2(n,l)]
  
  # Remainder of positive m (25)
  for l in range(0,lmax+1):
    for m in range(l-1, -1, -1):
      for n in range(l, -l, -1):
        vals[l][zos1(m,l),zos2(n-1,l)] = np.sqrt(float(l*(l+1)-m*(m+1))/(l*(l+1)-n*(n-1))) * vals[l][zos1(m+1,l),zos2(n,l)
            ] + float(m+n)/np.sqrt(float(l*(l+1)-n*(n-1))) * s / (1+c) * vals[l][zos1(m,l),zos2(n,l)]

  # Fill in the negative m (27)
  for l in range(0,lmax+1):
    for m in range(-l,0):
      for n in range(-l,l+1):
        vals[l][zos1(m,l),zos2(n,l)] = (-1)**(m+n) * vals[l][zos1(-m,l),zos2(-n,l)]

  # Adjust?
  """
  for l in range(0,lmax+1):
    for m in range(-l,l+1):
      for n in range(-l,l+1):
        vals[l][zos(m,l),zos(l,l)] /= ((-1)**(l-m) * (1j)**(n-m))
  """

  return vals

if __name__ == "__main__":
  mv = 3
  t = np.pi/3

  d1 = W_comp_all(t, mv)
  d2 = [np.zeros((2*i+1,2*i+1), dtype = complex) for i in range(mv+1)] 

  for l in range(mv+1):
    for m in range(-l,l+1):
      for n in range(-l,l+1):
        d2[l][m,n] = P(t,l,m,n)

  for l in range(mv+1):
    print(d1[l])
    print(d2[l])

  """
  for l in range(mv+1):
    for m in range(-l,l+1):
      for n in range(-l,l+1):
        print(l,m,n)
        print(P(t, l, m, n), P_alt(t,l,m,n))
        if np.abs(P(t,l,m,n) - P_alt(t,l,m,n)) > 1.0e-15:
          raise Exception()
  """
