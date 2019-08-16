from neml.math import rotations, nemlmath
from neml.cp import harmonics

import numpy as np
import numpy.polynomial as npp
from math import factorial

import unittest

class TestHarmonics(unittest.TestCase):
  def setUp(self):
    self.A = rotations.Orientation(30.0, 60.0, 80.0, 
        angle_type = "degrees", convention = "bunge")

  def test_SO3_associated(self):
    xs = [0.01,0.5,0.99]

    for n in range(7):
      for i in range(-n, n+1):
        for j in range(-n, n+1):
          for x in xs:
            P1 = P(x, n, i, j)
            P2 = harmonics.P_SO3(n, i, j, x)
            self.assertAlmostEqual(P2, P1)

class TestMetrics(unittest.TestCase):
  def setUp(self):
    self.A = rotations.Orientation(30.0, 60.0, 80.0, 
        angle_type = "degrees", convention = "bunge")

  def test_distance(self):
    self.assertTrue(np.isclose(self.A.distance(self.A),0.0))

  def test_ds(self):
    alpha = np.pi/3*2.5
    beta = np.pi/5*1.1
    gamma = np.pi*0.64

    q = rotations.Orientation(alpha, beta, gamma, convention = "hopf")

    self.assertAlmostEqual(
        np.sin(beta), np.sin(q.to_hopf()[1]))
    self.assertAlmostEqual(
        np.sin(beta), harmonics.ds(q))

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

def P(x, l, m, n):
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
