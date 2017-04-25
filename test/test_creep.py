import sys
sys.path.append('..')

from neml import creep
import unittest

from common import *

import numpy as np
import numpy.linalg as la

class CommonScalarCreep(object):
  """
    Common tests to impose on scalar creep laws.
  """
  def test_dg_ds(self):
    dfn = lambda x: self.model.g(x, self.e, self.t, self.T)
    nderiv = differentiate(dfn, self.s)
    cderiv = self.model.dg_ds(self.s, self.e, self.t, self.T)
    self.assertTrue(np.isclose(nderiv, cderiv, rtol = 1.0e-4))

  def test_dg_de(self):
    dfn = lambda x: self.model.g(self.s, x, self.t, self.T)
    nderiv = differentiate(dfn, self.e)
    cderiv = self.model.dg_de(self.s, self.e, self.t, self.T)
    self.assertTrue(np.isclose(nderiv, cderiv, rtol = 1.0e-4))

  def test_dg_dt(self):
    dfn = lambda x: self.model.g(self.s, self.e, x, self.T)
    nderiv = differentiate(dfn, self.t)
    cderiv = self.model.dg_dt(self.s, self.e, self.t, self.T)
    self.assertTrue(np.isclose(nderiv, cderiv))

  def test_dg_dT(self):
    dfn = lambda x: self.model.g(self.s, self.e, self.t, x)
    nderiv = differentiate(dfn, self.T)
    cderiv = self.model.dg_dT(self.s, self.e, self.t, self.T)
    self.assertTrue(np.isclose(nderiv, cderiv))

class TestPowerLawCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.A = 1.0e-6
    self.n = 5.0

    self.model = creep.PowerLawCreep(self.A, self.n) 

    self.T = 300.0
    self.e = 0.1
    self.s = 150.0
    self.t = 10.0

  def test_properties(self):
    self.assertTrue(np.isclose(self.A, self.model.A(self.T)))
    self.assertTrue(np.isclose(self.n, self.model.n(self.T)))

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    g_calc = self.A * self.s ** (self.n)
    self.assertTrue(np.isclose(g_direct, g_calc))

class TestNortonBaileyCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.A = 1.0e-6
    self.m = 0.25
    self.n = 5.0

    self.model = creep.NortonBaileyCreep(self.A, self.m, self.n) 

    self.T = 300.0
    self.e = 0.1
    self.s = 150.0
    self.t = 10.0

  def test_properties(self):
    self.assertTrue(np.isclose(self.A, self.model.A(self.T)))
    self.assertTrue(np.isclose(self.m, self.model.m(self.T)))
    self.assertTrue(np.isclose(self.n, self.model.n(self.T)))

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    g_calc = self.m * self.A**(1.0/self.m) * self.s**(self.n/self.m) * self.e**((self.m-1.0)/self.m)
    self.assertTrue(np.isclose(g_direct, g_calc))

