import sys
sys.path.append('..')

from neml import interpolate
import unittest

import numpy as np
import numpy.random as ra

class TestPolynomialInterpolate(unittest.TestCase):
  def setUp(self):
    self.n = 5
    self.coefs = list((ra.random((self.n,)) * 0.5 - 1.0) * 2.0)
    self.x = ra.random((1,))[0]
    self.interpolate = interpolate.PolynomialInterpolate(self.coefs)

  def test_interpolate(self):
    self.assertTrue(np.isclose(np.polyval(self.coefs, self.x),
      self.interpolate(self.x)))

class TestConstantInterpolate(unittest.TestCase):
  def setUp(self):
    self.v = ra.random((1,))[0]
    self.interpolate = interpolate.ConstantInterpolate(self.v)
    self.x = ra.random((1,))[0]

  def test_interpolate(self):
    self.assertTrue(np.isclose(self.v, self.interpolate(self.x)))

class TestMTSInterpolate(unittest.TestCase):
  def setUp(self):
    self.y0 = 100.0
    self.D = 50.0
    self.x0 = 200.0

    self.interpolate = interpolate.MTSShearInterpolate(self.y0, self.D, self.x0)

    self.x = 300.0

  def test_interpolate(self):
    should = self.y0 - self.D / (np.exp(self.x0 / self.x) - 1.0)
    act = self.interpolate(self.x)
    self.assertTrue(np.isclose(should, act))
