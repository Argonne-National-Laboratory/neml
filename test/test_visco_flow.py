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
    if (s > 0.0):
      v = s**self.n
    else:
      v = 0.0
    self.assertTrue(np.isclose(self.model.g(s), v))


class CommonFlowRule(object):
  """
    Tests common to all flow rules.
  """
  def gen_stress(self):
    s = np.array([200,0,100,0,-50,0])
    return s

  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_y(self):
    stress = self.gen_stress()
    hist = self.gen_hist()
    
    q = self.hrule.q(hist, self.T)
    fv = self.surface.f(stress, q, self.T)
    if fv > 0.0:
      should = self.g.g(fv) / self.eta
    else:
      should = 0.0
    
    ffm = self.model.y(stress, hist, self.T)
    
    self.assertTrue(np.isclose(should, ffm))


  def test_dy_ds(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.y(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dy_ds(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3, atol = 1.0e-3))

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

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

class TestPerzynaIsoJ2Voce(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.s0 = 180.0
    self.R = 150.0
    self.d = 10.0

    self.n = 5.0
    self.eta = 20.0

    self.surface = surfaces.IsoJ2()
    self.hrule = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)
    self.g = visco_flow.GPowerLaw(self.n)

    self.model = visco_flow.PerzynaFlowRule(self.surface, self.hrule, self.g, self.eta)

    self.hist0 = np.zeros((1,))
    self.T = 300.0

  def gen_hist(self):
    return np.array([0.01])

  def test_properties(self):
    self.assertTrue(np.isclose(self.eta, self.model.eta))


class TestPerzynaIsoJ2Linear(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.s0 = 100.0
    self.Kp = 1500.0

    self.n = 2.0
    self.eta = 100.0

    self.surface = surfaces.IsoJ2()
    self.hrule = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    self.g = visco_flow.GPowerLaw(self.n)

    self.model = visco_flow.PerzynaFlowRule(self.surface, self.hrule, self.g, self.eta)

    self.hist0 = np.zeros((1,))
    self.T = 300.0

  def gen_hist(self):
    return np.array([0.01])

  def test_properties(self):
    self.assertTrue(np.isclose(self.eta, self.model.eta))

class CommonFluidity(object):
  """
    Common fluidity model tests
  """
  def gen_hist(self):
    return ra.random((1,))[0]

  def test_deta(self):
    a = self.gen_hist()
    fm = self.model.deta(a)
    dfn = lambda a: self.model.eta(a)
    fn = differentiate(dfn, a)
    self.assertTrue(np.isclose(fm, fn))

class TestConstantFluidity(unittest.TestCase, CommonFluidity):
  def setUp(self):
    self.eta = 200.0
    self.model = visco_flow.ConstantFluidity(self.eta)

  def test_eta(self):
    a = self.gen_hist()
    self.assertTrue(np.isclose(self.eta, self.model.eta(a)))
