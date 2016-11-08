import sys
sys.path.append('..')

from neml import ri_flow, surfaces, hardening
from common import *

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

class CommonFlowRule(object):
  """
    Tests common to all flow rules.
  """
  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_df_ds(self):
    dfn = lambda s: self.model.f(s, self.hist_test, self.T_test)
    num = differentiate(dfn, self.stress_test)
    exact = self.model.df_ds(self.stress_test, self.hist_test, self.T_test)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_df_da(self):
    dfn = lambda a: self.model.f(self.stress_test, a, self.T_test)
    num = differentiate(dfn, self.hist_test)
    exact = self.model.df_da(self.stress_test, self.hist_test, self.T_test)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_ds(self):
    dfn = lambda s: self.model.g(s, self.hist_test, self.T_test)
    num = differentiate(dfn, self.stress_test)
    exact = self.model.dg_ds(self.stress_test, self.hist_test, self.T_test)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_da(self):
    dfn = lambda a: self.model.g(self.stress_test, a, self.T_test)
    num = differentiate(dfn, self.hist_test)
    exact = self.model.dg_da(self.stress_test, self.hist_test, self.T_test)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_ds(self):
    dfn = lambda s: self.model.h(s, self.hist_test, self.T_test)
    num = differentiate(dfn, self.stress_test)
    exact = self.model.dh_ds(self.stress_test, self.hist_test, self.T_test)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_da(self):
    dfn = lambda a: self.model.h(self.stress_test, a, self.T_test)
    num = differentiate(dfn, self.hist_test)
    exact = self.model.dh_da(self.stress_test, self.hist_test, self.T_test)
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))


class TestRateIndependentAssociativeFlowJ2Linear(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.s0 = 200.0
    self.K = 15000.0

    self.surface = surfaces.IsoJ2()
    self.hardening = hardening.LinearIsotropicHardeningRule(self.s0, self.K)
    self.model = ri_flow.RateIndependentAssociativeFlow(self.surface, self.hardening)

    self.hist0 = np.zeros((1,))
    
    self.stress_test = 100.0*ra.random((6,))
    self.hist_test = np.abs(ra.random((1,)))
    self.T_test = 300.0
  
  def test_flow_rule(self):
    dev_stress = self.stress_test - np.array([1,1,1,0,0,0]) * np.sum(self.stress_test[:3]) / 3.0
    should_be = dev_stress / la.norm(dev_stress)
    
    self.assertTrue(np.allclose(should_be, 
      self.model.g(self.stress_test, self.hist_test, self.T_test)))

  def test_hardening_rule(self):
    should_be = np.array([np.sqrt(2.0/3.0)])
    self.assertTrue(np.allclose(should_be,
      self.model.h(self.stress_test, self.hist_test, self.T_test)))
