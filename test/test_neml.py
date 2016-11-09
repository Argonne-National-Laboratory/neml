import sys
sys.path.append('..')

from neml import neml, elasticity, ri_flow, hardening, surfaces
from common import *

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

class CommonMatModel(object):
  """
    Tests that could be applied to all material models
  """
  def tests_tangent_proportional_strain(self):
    t_n = 0.0
    strain_n = np.zeros((6,))
    stress_n = np.zeros((6,))
    hist_n = np.copy(self.hist0)

    for i,m in enumerate(np.linspace(0,1,self.nsteps)):
      t_np1 = self.tfinal * m
      strain_np1 = self.efinal * m

      stress_np1, hist_np1, A_np1 = self.model.update_sd(strain_np1,
          strain_n, self.T, self.T, t_np1, t_n, stress_n, hist_n)
      dfn = lambda e: self.model.update_sd(e,
          strain_n, self.T, self.T, t_np1, t_n, stress_n, hist_n)[0]
      num_A = differentiate(dfn, strain_np1)
      
      self.assertTrue(np.allclose(num_A, A_np1, rtol = 1.0e-2))
      
      strain_n = strain_np1
      stress_n = stress_np1
      hist_n = hist_np1

class TestLinearElastic(CommonMatModel, unittest.TestCase):
  """
    Linear elasticity, as a benchmark
  """
  def setUp(self):
    self.hist0 = np.array([])
    self.E = 92000.0
    self.nu = 0.3

    self.mu = self.E/(2*(1+self.nu))
    self.K = self.E/(3*(1-2*self.nu))

    shear = elasticity.ConstantShearModulus(self.mu)
    bulk = elasticity.ConstantBulkModulus(self.K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)
    self.model = neml.SmallStrainElasticity(elastic)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

class TestRIAPlasticityJ2Linear(unittest.TestCase, CommonMatModel):
  """
    Test the rate-independent plasticity algorithm with a linearly
    isotropically hardening yield surface.
  """
  def setUp(self):
    self.hist0 = np.zeros((7,))

    self.E = 92000.0
    self.nu = 0.3

    self.mu = self.E/(2*(1+self.nu))
    self.K = self.E/(3*(1-2*self.nu))

    self.s0 = 180.0
    self.Kp = self.E / 100

    shear = elasticity.ConstantShearModulus(self.mu)
    bulk = elasticity.ConstantBulkModulus(self.K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoJ2()
    hrule = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model = neml.SmallStrainRateIndependentPlasticity(elastic, flow)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

  def test_residual(self):
    e_np1 = ra.random((6,))
    h_n = ra.random((7,))
    T = 300.0

    self.model.set_trial_state(e_np1, h_n, T)

    x = ra.random((8,))
    R, J = self.model.RJ(x)

    dfn = lambda y: self.model.RJ(y)[0]
    nJ = differentiate(dfn, x)

    self.assertTrue(np.allclose(J, nJ, rtol = 1.0e-3))
