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
      
      if not np.allclose(num_A, A_np1, rtol = 1e-3, atol = 1e-1):
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
    self.Kp = self.E/10

    self.elastic = elasticity.IsotropicLinearElasticModel(self.mu,
        "shear", self.K, "bulk")

    surface = surfaces.IsoJ2()
    hrule = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model = models.SmallStrainRateIndependentPlasticity(self.elastic,
        flow, max_divide = 5, force_divide = True)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

class TestDirectIntegrateChaboche(unittest.TestCase, CommonMatModel):
  """
    Test Chaboche's VP model with our new direct integrator
  """
  def setUp(self):
    n = 20.0
    eta = 108.0
    sY = 89.0

    Q = 165.0
    b = 12.0
    
    self.m = 3

    C1 = 80.0e3
    C2 = 14.02e3
    C3 = 3.333e3

    y1 = 0.9e3
    y2 = 1.5e3
    y3 = 1.0

    surface = surfaces.IsoKinJ2()
    iso = hardening.VoceIsotropicHardeningRule(sY, Q, b)
    cs = [C1, C2, C3]
    gs = [y1, y2, y3]
    As = [0.0, 0.0, 0.0]
    ns = [1.0, 1.0, 1.0]
    gmodels = [hardening.ConstantGamma(g) for g in gs]
    hmodel = hardening.Chaboche(iso, cs, gmodels, As, ns)

    fluidity = visco_flow.ConstantFluidity(eta)

    self.hist0 = np.zeros((19,))
    self.T = 300.0

    vmodel = visco_flow.ChabocheFlowRule(surface, hmodel, fluidity, n)

    E = 92000.0
    nu = 0.3

    mu = E/(2*(1+nu))
    K = E/(3*(1-2*nu))

    self.elastic = elasticity.IsotropicLinearElasticModel(mu,
        "shear", K, "bulk")

    flow = general_flow.TVPFlowRule(self.elastic, vmodel)

    self.model = models.GeneralIntegrator(self.elastic, flow,
        max_divide = 3, force_divide = True)

    self.efinal = np.array([0.05,0,0,0.02,0,-0.01])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 100
