import sys
sys.path.append('..')

from neml import interpolate, solvers, neml, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonMatModel(object):
  """
    Tests that could be applied to all material models
  """
  def tests_tangent_proportional_strain(self):
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
      
      if i != 0:
        self.assertTrue(np.allclose(num_A, A_np1, rtol = 1.0e-3, atol = 1.0e2))
      
      strain_n = strain_np1
      stress_n = stress_np1
      hist_n = hist_np1
      u_n = u_np1
      p_n = p_np1

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

    shear = elasticity.ShearModulus(self.mu)
    bulk = elasticity.BulkModulus(self.K)
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
    return np.array([0.01,-0.01,0.02,0.03,-0.005,0.01])

  def gen_T(self):
    return 500.0

  def gen_time(self):
    return 1.25

  def test_jacobian(self):
    e_np1 = self.gen_strain()
    T_np1 = self.gen_T()
    h_n = self.gen_hist()
    t_np1 = self.gen_time()

    ts = self.model.make_trial_state(e_np1, np.zeros((6,)), T_np1, T_np1,
        t_np1, 0.0, np.zeros((6,)), h_n)

    x = self.gen_x()

    R, J = self.model.RJ(x, ts)
    
    dfn = lambda y: self.model.RJ(y, ts)[0]
    nJ = differentiate(dfn, x)

    self.assertTrue(np.allclose(J, nJ, rtol = 1.0e-3))

    x = self.gen_x()

    R, J = self.model.RJ(x, ts)
    
    dfn = lambda y: self.model.RJ(y, ts)[0]
    nJ = differentiate(dfn, x)

    self.assertTrue(np.allclose(J, nJ, rtol = 1.0e-3))

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

    shear = elasticity.ShearModulus(self.mu)
    bulk = elasticity.BulkModulus(self.K)
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
    h = np.array(range(1,14)) / 15.0
    h[:6] = make_dev(h[:6])
    h[7:] = make_dev(h[7:])
    return h

  def gen_x(self):
    return np.array(list(self.gen_hist()) + list([0.25]))

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

    shear = elasticity.ShearModulus(self.mu)
    bulk = elasticity.BulkModulus(self.K)
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
    return np.array(range(1,8)) / 8.0

  def gen_x(self):
    return np.array(range(1,9)) / 9.0


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

    shear = elasticity.ShearModulus(self.mu)
    bulk = elasticity.BulkModulus(self.K)
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
    return np.array(range(1,8)) / 8.0

  def gen_x(self):
    return np.array(range(1,9)) / 9.0


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

    shear = elasticity.ShearModulus(self.mu)
    bulk = elasticity.BulkModulus(self.K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    hmodel = hardening.Chaboche(iso, self.cs, self.rs)

    flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hmodel)

    self.model = neml.SmallStrainRateIndependentPlasticity(elastic, flow,
        check_kt = False)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 30

  def gen_hist(self):
    h = np.array(range(1,6+1+6*self.n+1)) / (8*self.n)
    h[:6] = make_dev(h[:6])
    for i in range(self.n):
      h[7+i*6:7+(i+1)*6] = make_dev(h[7+i*6:7+(i+1)*6])
    return h

  def gen_x(self):
    return np.array(list(self.gen_hist()) + [0.1])

class TestDirectIntegrateCheboche(unittest.TestCase, CommonMatModel, CommonJacobian):
  """
    Test Cheboche's VP model with our new direct integrator
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
    cs = np.array([C1, C2, C3])
    gs = np.array([y1, y2, y3])
    hmodel = hardening.Chaboche(iso, cs, gs)

    fluidity = visco_flow.ConstantFluidity(eta)

    self.hist0 = np.zeros((19,))
    self.T = 300.0

    vmodel = visco_flow.ChabocheFlowRule(surface, hmodel, fluidity, n)

    E = 92000.0
    nu = 0.3

    mu = E/(2*(1+nu))
    K = E/(3*(1-2*nu))

    shear = elasticity.ShearModulus(mu)
    bulk = elasticity.BulkModulus(K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    flow = general_flow.TVPFlowRule(elastic, vmodel)

    self.model = neml.GeneralIntegrator(flow)

    self.efinal = np.array([0.05,0,0,0.02,0,-0.01])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 100

  def gen_hist(self):
    h = np.array([0.05,20,-30,40.0,5.0,2.0,40.0,-10,-20,10,50,30,10,-50,60,70,15,-15,10])
    h[1:7] = make_dev(h[1:7])
    h[7:13] = make_dev(h[7:13])
    h[13:19] = make_dev(h[13:19])
    return h

  def gen_x(self):
    x = [100.0,150.0,-300.0,-10.0,50.0,100.0] + list(self.gen_hist()*1.1)
    return np.array(x)

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

    shear = elasticity.ShearModulus(self.mu)
    bulk = elasticity.BulkModulus(self.K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoJ2()
    hrule = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)
    g = visco_flow.GPowerLaw(self.n)

    vmodel = visco_flow.PerzynaFlowRule(surface, hrule, g, self.eta)

    flow = general_flow.TVPFlowRule(elastic, vmodel)
    self.model = neml.GeneralIntegrator(flow)

    self.efinal = np.array([0.05,0,0,0.02,0,-0.01])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 100 

  def gen_hist(self):
    return np.array(range(1,8)) / 8.0

  def gen_x(self):
    return np.array(range(1,8)) / 8.0

class TestYaguchi(unittest.TestCase, CommonMatModel, CommonJacobian):
  """
    Test Cheboche's VP model with our new direct integrator
  """
  def setUp(self):
    vmodel = visco_flow.YaguchiGr91FlowRule()

    E_poly = np.array([-3.077e-1, 3.003e2, 1.269e5])
    nu = 0.3

    mu_poly = list(E_poly/(2*(1+nu)))
    K_poly = list(E_poly/(3*(1-2*nu)))

    shear = elasticity.ShearModulus(interpolate.PolynomialInterpolate(mu_poly))
    bulk = elasticity.BulkModulus(interpolate.PolynomialInterpolate(K_poly))
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    flow = general_flow.TVPFlowRule(elastic, vmodel)

    self.model = neml.GeneralIntegrator(flow)

    self.efinal = np.array([0.05,0,0,0.02,0,-0.01])
    self.tfinal = 10.0
    self.T = 550.0
    self.nsteps = 100

  def gen_hist(self):
    h = np.array([20,-30,40.0,5.0,2.0,40.0,-10,-20,10,50,30,10,10.0,5.0])
    h[:6] = make_dev(h[:6])
    h[6:12] = make_dev(h[6:12])
    return h

  def gen_x(self):
    x = [100.0,150.0,-300.0,-10.0,50.0,100.0] + list(self.gen_hist()*1.1)
    return np.array(x)
