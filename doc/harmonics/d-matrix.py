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

if __name__ == "__main__":
  mv = 3
  t = np.pi/3

  for l in range(mv+1):
    for m in range(-l,l+1):
      for n in range(-l,l+1):
        print(l,m,n)
        print(P(t, l, m, n), P_alt(t,l,m,n))
        if np.abs(P(t,l,m,n) - P_alt(t,l,m,n)) > 1.0e-15:
          raise Exception()
