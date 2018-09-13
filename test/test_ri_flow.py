import sys
sys.path.append('..')

from neml import ri_flow, surfaces, hardening
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonFlowRule(object):
  """
    Tests common to all flow rules.
  """
  def gen_stress(self):
    s = np.array([50.0,75.0,-50.0,25.0,75.0,-25.0])
    return s

  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_df_ds(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.f(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.df_ds(stress, hist, self.T)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_df_da(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.f(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.df_da(stress, hist, self.T)
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


class TestRateIndependentAssociativeFlowJ2Linear(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.s0 = 200.0
    self.K = 15000.0

    self.surface = surfaces.IsoJ2()
    self.hardening = hardening.LinearIsotropicHardeningRule(self.s0, self.K)
    self.model = ri_flow.RateIndependentAssociativeFlow(self.surface, self.hardening)

    self.hist0 = np.zeros((1,))
    
    self.T = 300.0

  def gen_hist(self):
    return 0.25
  
  def test_flow_rule(self):
    stress = self.gen_stress()
    dev_stress = make_dev(stress)
    should_be = dev_stress / la.norm(dev_stress)
    
    self.assertTrue(np.allclose(should_be, 
      self.model.g(stress, self.gen_hist(), self.T)))

  def test_hardening_rule(self):
    should_be = np.array([np.sqrt(2.0/3.0)])
    self.assertTrue(np.allclose(should_be,
      self.model.h(self.gen_stress(), self.gen_hist(), self.T)))


class TestRIChabocheLinear(unittest.TestCase, CommonFlowRule):
  """
    Test the rate independent Chaboche model with linear isotropic hardening
  """
  def setUp(self):
    self.s0 = 200.0
    self.K = 15000.0

    self.n = 4
    self.cs = range(1,self.n+1)
    self.rs = range(1, self.n+1)

    self.iso = hardening.LinearIsotropicHardeningRule(self.s0, self.K)
    self.gmodels = [hardening.ConstantGamma(g) for g in self.rs]
    As = [0.0] * self.n
    ns = [1.0] * self.n
    self.hardening = hardening.Chaboche(self.iso, self.cs, self.gmodels, As, 
        ns)
    self.surface = surfaces.IsoKinJ2()

    self.model = ri_flow.RateIndependentNonAssociativeHardening(self.surface, self.hardening)

    self.hist0 = np.zeros((1+6*self.n,))

    self.T = 300.0

  def gen_hist(self):
    hist = np.array(range(1,2+self.n*6))  / (self.n*7)
    hist[1:] = (1.0 - 2.0 * hist[1:]) * 100.0
    for i in range(self.n):
      hist[1+i*6:1+(i+1)*6] = make_dev(hist[1+i*6:1+(i+1)*6])
    return hist

  def test_flow_rule(self):
    stress = self.gen_stress()
    s = make_dev(stress)
    hist = self.gen_hist()

    X = sum(hist[1+i*6:1+(i+1)*6] for i in range(self.n))

    n = (s+X) / la.norm(s+X)

    self.assertTrue(np.allclose(n,
      self.model.g(stress, hist, self.T)))
