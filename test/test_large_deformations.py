import sys
sys.path.append('..')

from neml import models, elasticity, surfaces, hardening, ri_flow
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonLD(object):
  """
    Generic tests that can be applied to arbitrary models
  """
  def test_tangent_at_points(self):
    for c in self.conditions:
      hist_n = c['hist_n']
      stress_n = c['stress_n']
      d_n = c['d_n']
      d_np1 = c['d_np1']
      w_n = c['w_n']
      w_np1 = c['w_np1']
      t_n = 0.0
      t_np1 = c['dt']
      T_n = c['T']
      T_np1 = c['T']
      u_n = 0.0
      p_n = 0.0

      stress_np1, hist_np1, A_np1, B_np1, u_np1, p_np1 = self.model.update_ld_inc(
          d_np1, d_n, w_np1, w_n, T_np1, T_n, t_np1, t_n, stress_n, hist_n,
          u_n, p_n)

      ufn_d = lambda d: self.model.update_ld_inc(
          d, d_n, w_np1, w_n, T_np1, T_n, t_np1, t_n, stress_n, hist_n,
          u_n, p_n)[0]
      ufn_w = lambda w: self.model.update_ld_inc(
          d_np1, d_n, w, w_n, T_np1, T_n, t_np1, t_n, stress_n, hist_n,
          u_n, p_n)[0]

      A_num = differentiate(ufn_d, d_np1)
      B_num = differentiate(ufn_w, w_np1)
      
      print("A")
      print((A_np1-A_num)/A_num)

      print("B")
      print((B_np1-B_num)/B_num)

      self.assertTrue(np.allclose(A_np1, A_num, rtol = 1.0e-2))
      self.assertTrue(np.allclose(B_np1, B_num, rtol = 1.0e-2))

  def test_tangent_proportional(self):
    for dirs in self.directions:
      ddir = dirs['d']
      wdir = dirs['w']

      t_n = 0.0
      hist_n = self.model.init_store()
      d_n = np.zeros((6,))
      w_n = np.zeros((3,))
      stress_n = np.zeros((6,))

      for i in range(self.nsteps):
        t_np1 = t_n + self.dt
        d_np1 = d_n + ddir
        w_np1 = w_n + wdir

        stress_np1, hist_np1, A_np1, B_np1, u_np1, p_np1 = self.model.update_ld_inc(
            d_np1, d_n, w_np1, w_n, self.T, self.T, t_np1, t_n, stress_n, hist_n,
            0.0, 0.0)

        ufn_d = lambda d: self.model.update_ld_inc(
            d, d_n, w_np1, w_n, self.T, self.T, t_np1, t_n, stress_n, hist_n,
            0.0, 0.0)[0]
        ufn_w = lambda w: self.model.update_ld_inc(
            d_np1, d_n, w, w_n, self.T, self.T, t_np1, t_n, stress_n, hist_n,
            0.0, 0.0)[0]

        A_num = differentiate(ufn_d, d_np1)
        B_num = differentiate(ufn_w, w_np1)

        self.assertTrue(np.allclose(A_num, A_np1, rtol = 1.0e-2))
        self.assertTrue(np.allclose(B_num, B_np1, rtol = 1.0e-2))

        t_n = t_np1
        d_n = np.copy(d_np1)
        w_n = np.copy(w_np1)
        stress_n = np.copy(stress_np1)
        hist_n = np.copy(hist_np1)


class TestLinearElastic(CommonLD, unittest.TestCase):
  def setUp(self):
    self.E = 100000.0
    self.nu = 0.29

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E,
        "youngs", self.nu, "poissons")
    self.model = models.SmallStrainElasticity(self.elastic)
    
    self.conditions = [
        {'hist_n': np.zeros((6,)),
          'stress_n': np.zeros((6,)),
          'd_n': np.zeros((6,)),
          'd_np1': np.array([0.1,0.05,-0.025,0.15,0.2,-0.05]),
          'w_n': np.zeros((3,)),
          'w_np1': np.array([-0.15,0.1,0.05]),
          'dt': 1.0,
          'T': 300.0},
        {'hist_n': np.zeros((6,)),
          'stress_n': np.array([10.0,-5.0,30.0,-5.0,10.0,15.0]),
          'd_n': np.zeros((6,)),
          'd_np1': np.array([0.05,0.5,0.25,0.20,-0.2,0.25]),
          'w_n': np.zeros((3,)),
          'w_np1': np.array([-0.15,0.25,0.05]),
          'dt': 1.0,
          'T': 300.0}
        ]

    self.directions = [
        {'d': np.array([0.1,-0.15,0.2,-0.05,0.15,0.25]),
          'w': np.array([0.25,-0.15,0.1])}]

    self.nsteps = 10
    self.dt = 1.0
    self.T = 300.0

class TestSimplePlastic(CommonLD, unittest.TestCase):
  def setUp(self):
    self.E = 100000.0
    self.nu = 0.29
    self.sY = 100.0
    self.H = 1200.0

    elastic = elasticity.IsotropicLinearElasticModel(self.E,
        "youngs", self.nu, "poissons")
    surf = surfaces.IsoJ2()
    hard = hardening.LinearIsotropicHardeningRule(self.sY, self.H)
    flow = ri_flow.RateIndependentAssociativeFlow(surf, hard)
    self.model = models.SmallStrainRateIndependentPlasticity(elastic, flow)

    self.conditions = [
        {'hist_n': self.model.init_store(),
          'stress_n': np.zeros((6,)),
          'd_n': np.zeros((6,)),
          'd_np1': np.array([0.1,0.05,-0.025,0.15,0.2,-0.05]),
          'w_n': np.zeros((3,)),
          'w_np1': np.array([-0.15,0.1,0.05]),
          'dt': 1.0,
          'T': 300.0},
        {'hist_n': self.model.init_store(),
          'stress_n': np.array([10.0,-5.0,30.0,-5.0,10.0,15.0]),
          'd_n': np.zeros((6,)),
          'd_np1': np.array([0.05,0.5,0.25,0.20,-0.2,0.25]),
          'w_n': np.zeros((3,)),
          'w_np1': np.array([-0.15,0.25,0.05]),
          'dt': 1.0,
          'T': 300.0}
        ]

    self.directions = [
        {'d': np.array([0.1,-0.15,0.2,-0.05,0.15,0.25]),
          'w': np.array([0.25,-0.15,0.1])}]

    self.nsteps = 10
    self.dt = 1.0
    self.T = 300.0
