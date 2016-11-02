import sys
sys.path.append('..')

from neml.hardening import *
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonAssociative(object):
  """
    Inherit from this to test associative hardening rules.
  """
  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.alpha0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.alpha0))

  def test_D(self):
    dfn = lambda a: self.model.q(a, self.T)
    n_d = differentiate(dfn, self.alpha)

    self.assertTrue(np.allclose(self.model.D(self.alpha, self.T), n_d))

  def test_D_inv(self):
    D = self.model.D(self.alpha, self.T)
    self.assertTrue(np.allclose(la.inv(D), self.model.D_inv(self.alpha, self.T)))

class TestIsoJ2LinearAHardening(unittest.TestCase, CommonAssociative):
  def setUp(self):
    self.alpha = np.array([0.1])
    self.alpha0 = np.zeros((1,))

    self.K0 = 100.0
    self.Kp = 1001.0

    self.T = 300.0

    self.model = IsoJ2LinearAHardening(self.K0, self.Kp)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.K0, self.K0))
    self.assertTrue(np.isclose(self.model.Kp, self.Kp))

  def test_q(self):
    self.assertTrue(np.allclose(self.model.q(self.alpha, self.T), -self.K0 - self.Kp * self.alpha[0]))
