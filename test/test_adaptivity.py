import sys
sys.path.append('..')

from neml import interpolate, solvers, models, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonMatModel(object):
  """
    Tests that could be applied to all material models
  """
  def test_tangent_proportional_strain(self):
    t_n = 0.0
    strain_n = np.zeros((6,))
    stress_n = np.zeros((6,))
    hist_n = self.model.init_store()

    u_n = 0.0
    p_n = 0.0

    for i,m in enumerate(np.linspace(0,1,self.nsteps)):
      t_np1 = self.tfinal * m
      strain_np1 = self.efinal * m
      
      stress_np1, hist_np1, A_np1, u_np1, p_np1 = self.model.update_sd(
          strain_np1, strain_n, self.T, self.T, t_np1, t_n, stress_n, hist_n,
          u_n, p_n)
      dfn = lambda e: self.model.update_sd(e,
          strain_n, self.T, self.T, t_np1, t_n, stress_n, hist_n, u_n, p_n)[0]
      num_A = differentiate(dfn, strain_np1, eps = 1.0e-9)

      print(A_np1)
      print(num_A)
      
      self.assertTrue(np.allclose(num_A, A_np1, rtol = 1.0e-3, atol = 1.0e-1))
      
      strain_n = strain_np1
      stress_n = stress_np1
      hist_n = hist_np1
      u_n = u_np1
      p_n = p_np1

class TestPerfectPlasticity(unittest.TestCase, CommonMatModel):
  def setUp(self):
    self.hist0 = np.zeros((0,))

    self.E = 92000.0
    self.nu = 0.3

    self.mu = self.E/(2*(1+self.nu))
    self.K = self.E/(3*(1-2*self.nu))

    self.s0 = 180.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.mu,
        "shear", self.K, "bulk")

    surface = surfaces.IsoJ2()
    self.model = models.SmallStrainPerfectPlasticity(
        self.elastic, surface, self.s0,
        max_divide = 5, force_divide = True)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10
