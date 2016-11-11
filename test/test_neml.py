import sys
sys.path.append('..')

from neml import neml, elasticity, ri_flow, hardening, surfaces, visco_flow
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

class CommonJacobian(object):
  """
    Common jacobian tests
  """
  def gen_strain(self):
    return ra.random((6,))

  def gen_T(self):
    return ra.random((1,))[0]*1000.0

  def gen_time(self):
    return ra.random((1,))[0] * 10.0

  def test_jacobian(self):
    e_np1 = self.gen_strain()
    T_np1 = self.gen_T()
    h_n = self.gen_hist()
    t_np1 = self.gen_time()

    self.model.set_trial_state(e_np1, h_n, T_np1, t_np1, 0.0)

    x = self.gen_x()

    R, J = self.model.RJ(x)
    
    dfn = lambda y: self.model.RJ(y)[0]
    nJ = differentiate(dfn, x)
    
    self.assertTrue(np.allclose(J, nJ, rtol = 1.0e-2))

class TestRIAPlasticityCombinedLinearLinear(unittest.TestCase, CommonMatModel, CommonJacobian):
  """
    Test combined isotropic/kinematic hardening with linear rules
  """
  def setUp(self):
    self.hist0 = np.zeros((13,))

    self.E = 92000.0
    self.nu = 0.3

    self.mu = self.E/(2*(1+self.nu))
    self.K = self.E/(3*(1-2*self.nu))

    self.s0 = 180.0
    self.Kp = 1000.0

    self.H = 1000.0

    shear = elasticity.ConstantShearModulus(self.mu)
    bulk = elasticity.ConstantBulkModulus(self.K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model = neml.SmallStrainRateIndependentPlasticity(elastic, flow, check_kt = False)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

  def gen_hist(self):
    h = ra.random((13,))
    h[:6] = make_dev(h[:6])
    h[7:] = make_dev(h[7:])
    return h

  def gen_x(self):
    return np.array(list(self.gen_hist()) + list(ra.random((1,))))

class TestRIAPlasticityJ2Linear(unittest.TestCase, CommonMatModel, CommonJacobian):
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

    self.model = neml.SmallStrainRateIndependentPlasticity(elastic, flow, check_kt = False)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

  def gen_hist(self):
    return ra.random((7,))

  def gen_x(self):
    return ra.random((8,))


class TestRIAPlasticityJ2Voce(unittest.TestCase, CommonMatModel, CommonJacobian):
  """
    Test the rate-independent plasticity algorithm with a Voce
    isotropically hardening yield surface.
  """
  def setUp(self):
    self.hist0 = np.zeros((7,))

    self.E = 92000.0
    self.nu = 0.3

    self.mu = self.E/(2*(1+self.nu))
    self.K = self.E/(3*(1-2*self.nu))

    self.s0 = 180.0
    self.R = 150.0
    self.d = 10.0

    shear = elasticity.ConstantShearModulus(self.mu)
    bulk = elasticity.ConstantBulkModulus(self.K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoJ2()
    hrule = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model = neml.SmallStrainRateIndependentPlasticity(elastic, flow, check_kt = False)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

  def gen_hist(self):
    return ra.random((7,))

  def gen_x(self):
    return ra.random((8,))


class TestRIChebocheLinear(unittest.TestCase, CommonMatModel, CommonJacobian):
  """
    Test Cheboche with linear isotropic hardening
  """
  def setUp(self):
    self.hist0 = np.zeros((13,))

    self.E = 92000.0
    self.nu = 0.3

    self.mu = self.E/(2*(1+self.nu))
    self.K = self.E/(3*(1-2*self.nu))

    self.s0 = 180.0
    self.Kp = 1000.0

    self.n = 2
    self.cs = np.array([10.0, 2.0])
    self.rs = np.array([5.0, 1.0])

    shear = elasticity.ConstantShearModulus(self.mu)
    bulk = elasticity.ConstantBulkModulus(self.K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    hmodel = hardening.Chaboche(iso, self.cs, self.rs)

    flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hmodel)

    self.model = neml.SmallStrainRateIndependentPlasticity(elastic, flow, check_kt = False)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 30

  def gen_hist(self):
    h = ra.random((6+1+6*self.n,))
    h[:6] = make_dev(h[:6])
    for i in range(self.n):
      h[7+i*6:7+(i+1)*6] = make_dev(h[7+i*6:7+(i+1)*6])
    return h

  def gen_x(self):
    return np.array(list(self.gen_hist()) + list(ra.random((1,))))


class TestPerzynaJ2Voce(unittest.TestCase, CommonMatModel, CommonJacobian):
  """
    Perzyna associated viscoplasticity w/ voce kinematic hardening
  """
  def setUp(self):
    self.hist0 = np.zeros((7,))

    self.hist0 = np.zeros((7,))

    self.E = 92000.0
    self.nu = 0.3

    self.mu = self.E/(2*(1+self.nu))
    self.K = self.E/(3*(1-2*self.nu))

    self.s0 = 180.0
    self.R = 150.0
    self.d = 10.0

    self.n = 2.0
    self.eta = 200.0

    shear = elasticity.ConstantShearModulus(self.mu)
    bulk = elasticity.ConstantBulkModulus(self.K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoJ2()
    hrule = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)
    g = visco_flow.GPowerLaw(self.n)

    flow = visco_flow.PerzynaFlowRule(surface, hrule, g, self.eta)
    self.model = neml.SmallStrainViscoPlasticity(elastic, flow)

    self.efinal = np.array([0.01,-0.005,0.02,-0.03,0.01,-0.015])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 100 

  def gen_hist(self):
    return ra.random((7,))

  def gen_x(self):
    return ra.random((8,))
