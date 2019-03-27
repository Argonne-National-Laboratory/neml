import sys
sys.path.append('..')

from neml import models, elasticity
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
      print(A_np1)
      print(A_num)

      print("B")
      print(B_np1)
      print(B_num)

      self.assertTrue(np.allclose(A_np1, A_num, rtol = 1.0e-3))
      self.assertTrue(np.allclose(B_np1, B_num, rtol = 1.0e-3))


class TestLinearElastic(CommonLD, unittest.TestCase):
  def setUp(self):
    self.E = 100000.0
    self.nu = 0.29

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E,
        "youngs", self.nu, "poissons")
    self.model = models.SmallStrainElasticity(self.elastic)
    
    self.conditions = [
        {'hist_n': np.array([]),
          'stress_n': np.zeros((6,)),
          'd_n': np.zeros((6,)),
          'd_np1': np.array([0.1,0.05,-0.025,0.15,0.2,-0.05]),
          'w_n': np.zeros((3,)),
          'w_np1': np.array([-0.15,0.1,0.05]),
          'dt': 1.0,
          'T': 300.0}
        ]
