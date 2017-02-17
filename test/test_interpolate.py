import sys
sys.path.append('..')

from neml import interpolate
import unittest

import numpy as np
import numpy.random as ra
import scipy.interpolate as inter

class TestPolynomialInterpolate(unittest.TestCase):
  def setUp(self):
    self.n = 5
    self.coefs = list((ra.random((self.n,)) * 0.5 - 1.0) * 2.0)
    self.x = ra.random((1,))[0]
    self.interpolate = interpolate.PolynomialInterpolate(self.coefs)

  def test_interpolate(self):
    self.assertTrue(np.isclose(np.polyval(self.coefs, self.x),
      self.interpolate(self.x)))

class TestPiecewiseLinearInterpolate(unittest.TestCase):
  def setUp(self):
    self.validx = [-10.0, -2.0, 1.0, 2.0, 5.0, 15.0]
    self.invalidx_1 = [-2.0, -10.0, 1.0, 2.0, 5.0, 15.0]
    self.invalidx_2 = [-10.0, -2.0, 1.0, 2.0, 5.0]
    self.points = [50.0, -25.0, 1.0, 0.0, 0.0, 10.0]

    self.valid = interpolate.PiecewiseLinearInterpolate(self.validx, 
        self.points)
    self.invalid1 = interpolate.PiecewiseLinearInterpolate(self.invalidx_1, 
        self.points)
    self.invalid2 = interpolate.PiecewiseLinearInterpolate(self.invalidx_2, 
        self.points)

  def test_valid(self):
    self.assertTrue(self.valid.valid)
    self.assertFalse(self.invalid1.valid)
    self.assertFalse(self.invalid2.valid)

  def test_interpolate(self):
    testinter = inter.interp1d(self.validx, self.points,
        fill_value = (self.points[0], self.points[-1]),
        bounds_error = False)
    xs = np.linspace(-15.0,20.0)
    ys1 = [self.valid(x) for x in xs]
    ys2 = testinter(xs)
    self.assertTrue(np.allclose(ys1, ys2))

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
