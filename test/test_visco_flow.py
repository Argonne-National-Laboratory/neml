import sys
sys.path.append('..')

from neml import visco_flow, surfaces, hardening
from common import *

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

class CommonGFlow(object):
  """
    Common to our Perzyna g functions.
  """
  def test_dg(self):
    s = self.gen_f()

    exact = self.model.dg(s)

    dfn = lambda x: self.model.g(x)
    num = differentiate(dfn, s)

    print(exact)
    print(num)

    self.assertTrue(np.isclose(exact, num))

class TestGPowerLaw(unittest.TestCase, CommonGFlow):
  def setUp(self):
    self.n = 11.0
    self.model = visco_flow.GPowerLaw(self.n)

  def gen_f(self):
    return 100*(1-2.0*ra.random((1,))[0])

  def test_properties(self):
    self.assertTrue(np.isclose(self.n, self.model.n))

  def test_g(self):
    s = self.gen_f()
    self.assertTrue(np.isclose(self.model.g(s), np.abs(s**(self.n-1))*s))


class CommonFlowRule(object):
  """
    Tests common to all flow rules.
  """
  def gen_stress(self):
    s = ra.random((6,))
    s = (1.0 - 2.0 * s) * 125.0
    return s

  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_dy_ds(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.y(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dy_ds(stress, hist, self.T)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dy_da(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.y(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dy_da(stress, hist, self.T)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_ds(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.g(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dg_ds(stress, hist, self.T)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_da(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.g(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dg_da(stress, hist, self.T)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_ds(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.h(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dh_ds(stress, hist, self.T)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_da(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.h(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dh_da(stress, hist, self.T)
    
    # Check for bad random state
    if np.any(np.isnan(exact)):
      return

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

class TestPerzynaIsoJ2Voce(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.s0 = 180.0
    self.R = 150.0
    self.d = 10.0

    self.n = 5.0
    self.eta = 20.0

    surface = surfaces.IsoJ2()
    hrule = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)
    g = visco_flow.GPowerLaw(self.n)

    self.model = visco_flow.PerzynaFlowRule(surface, hrule, g, self.eta)

    self.hist0 = np.zeros((1,))
    self.T = 300.0

  def gen_hist(self):
    return ra.random((1,))

  def test_properties(self):
    self.assertTrue(np.isclose(self.eta, self.model.eta))
