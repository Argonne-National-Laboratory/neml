import sys
sys.path.append('..')

from neml import hardening
import unittest

from common import *

import numpy as np
import numpy.linalg as la
import numpy.random as ra


class CommonHardening(object):
  """
    Tests that can apply to all hardening rules
  """
  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_gradient(self):
    dfn = lambda x: self.model.q(x, self.T)
    ngrad = differentiate(dfn, self.hist_trial)
    self.assertTrue(np.allclose(ngrad, self.model.dq_da(self.hist_trial, self.T)))

class TestLinearIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0

    self.hist0 = np.array([0.0])
    
    self.hist_trial = np.abs(ra.random((1,)))
    self.T = 300.0

    self.model = hardening.LinearIsotropicHardeningRule(self.s0, self.K)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.s0, self.s0))
    self.assertTrue(np.isclose(self.model.K, self.K))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      np.array([-self.s0 - self.K * self.hist_trial[0]])))

class TestVoceIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.R = 100.0
    self.d = 10.0

    self.hist0 = np.array([0.0])
    
    self.hist_trial = np.abs(ra.random((1,)))
    self.T = 300.0

    self.model = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.s0, self.s0))
    self.assertTrue(np.isclose(self.model.R, self.R))
    self.assertTrue(np.isclose(self.model.d, self.d))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      np.array([-self.s0 - self.R * (1 - np.exp(-self.d*self.hist_trial[0]))])))

class TestLinearKinematicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.H = 1000.0

    self.hist0 = np.zeros((6,))
    
    self.hist_trial = ra.random((6,)) * 100
    self.hist_trial = self.hist_trial - np.array([1,1,1,0,0,0]) * sum(self.hist_trial[:3]) / 3.0
    self.T = 300.0

    self.model = hardening.LinearKinematicHardeningRule(self.H)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.H, self.H))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      -self.hist_trial * self.H))


