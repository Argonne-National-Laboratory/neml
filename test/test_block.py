import sys
sys.path.append('..')

from neml import models, block, elasticity, surfaces, hardening, ri_flow
from common import *

import unittest
import numpy as np
import numpy.random as ra

class TestBlockEvaluate(unittest.TestCase):
  def setUp(self):
    E = 150000.0
    nu = 0.3

    ys = 100.0
    H = 1000.0

    elastic = elasticity.IsotropicLinearElasticModel(E, "youngs",
        nu, "poissons")

    surface = surfaces.IsoJ2()
    hrule = hardening.LinearIsotropicHardeningRule(ys, H)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model = models.SmallStrainRateIndependentPlasticity(elastic, flow)
    
    self.nblock = 100

  def test_batch(self):
    e_np1 = np.zeros((self.nblock,3,3))
    e_np1[:,0,0] = np.linspace(0, 0.1, self.nblock)
    e_n = np.zeros((self.nblock,3,3))
    T_np1 = np.zeros((self.nblock,))
    T_n = np.zeros((self.nblock,))
    t_np1 = 1.0
    t_n = 0.0
    s_np1 = np.zeros((self.nblock,3,3))
    s_n = np.zeros((self.nblock,3,3))
    h_np1 = np.zeros((self.nblock, self.model.nstore))
    h_n = np.zeros((self.nblock, self.model.nstore))
    A_np1 = np.zeros((self.nblock,3,3,3,3))
    u_np1 = np.zeros((self.nblock,))
    u_n = np.zeros((self.nblock,))
    p_np1 = np.zeros((self.nblock,))
    p_n = np.zeros((self.nblock,))

    # Try with the block eval function
    block.block_evaluate(self.model,
        e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n, h_np1, h_n,
        A_np1, u_np1, u_n, p_np1, p_n)

    # Now compare to item-by-item evaluation
    for i in range(self.nblock):
      s, h, A, u, p = self.model.update_sd(
          sym(e_np1[i]), sym(e_n[i]), T_np1[i], T_n[i], 
          t_np1, t_n, sym(s_n[i]), h_n[i], u_n[i], 
          p_n[i])
      print(e_np1[i])
      print(usym(s))
      print(s_np1[i])
      self.assertTrue(np.allclose(usym(s), s_np1[i]))
      self.assertTrue(np.allclose(h, h_np1[i]))
      self.assertTrue(np.allclose(ms2ts(A), A_np1[i]))
      self.assertTrue(np.isclose(u_np1[i], u))
      self.assertTrue(np.isclose(p_np1[i], p))

class TestBlockConverters(unittest.TestCase):
  def setUp(self):
    self.n = 100

  def test_t2m(self):
    T = ra.random((self.n,3,3))
    T += T.transpose(0,2,1)

    A1 = np.zeros((self.n,6))
    for i in range(self.n):
      A1[i] = sym(T[i])

    A2 = block.t2m(T)

    self.assertTrue(np.allclose(A1,A2))

  def test_m2t(self):
    M = ra.random((self.n,6))

    A1 = np.zeros((self.n,3,3))
    for i in range(self.n):
      A1[i] = usym(M[i])
    A2 = block.m2t(M)

    self.assertTrue(np.allclose(A1,A2))

  def test_m42t4(self):
    M = ra.random((self.n,6,6))

    A1 = np.zeros((self.n,3,3,3,3))
    for i in range(self.n):
      A1[i] = ms2ts(M[i])

    A2 = block.m42t4(M)

    self.assertTrue(np.allclose(A1,A2))
