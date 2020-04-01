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

  def test_elastic_proportional_strain(self):
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
      
      estrain = self.model.elastic_strains(stress_np1, self.T, hist_np1)
      S = self.elastic.S(self.T)
      estrain_num = np.dot(S, stress_np1)
      self.assertTrue(np.allclose(estrain, estrain_num))
      
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

    self.elastic = elasticity.IsotropicLinearElasticModel(self.mu,
        "shear", self.K, "bulk")
    self.model = models.SmallStrainElasticity(self.elastic)

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

  def gen_start_strain(self):
    return np.zeros((6,))

  def gen_start_time(self):
    return 0.0

  def gen_start_stress(self):
    return np.zeros((6,))

  def test_jacobian(self):
    e_np1 = self.gen_strain()
    T_np1 = self.gen_T()
    h_n = self.gen_hist()
    t_np1 = self.gen_time()

    ts = self.model.make_trial_state(e_np1, self.gen_start_strain(), T_np1, T_np1,
        t_np1, self.gen_start_time(), self.gen_start_stress(), h_n)

    x = self.gen_x()

    R, J = self.model.RJ(x, ts)
    
    dfn = lambda y: self.model.RJ(y, ts)[0]
    nJ = differentiate(dfn, x)
   
    self.assertTrue(np.allclose(J, nJ, rtol = 1.0e-3, atol = 1e-1))

class TestPerfectPlasticity(unittest.TestCase, CommonMatModel, CommonJacobian):
  """
    Test J2 perfect plasticity
  """
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
    self.model = models.SmallStrainPerfectPlasticity(self.elastic, surface, 
        self.s0)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.ys(self.T), self.s0))

  def gen_hist(self):
    return np.zeros((0,))

  def gen_x(self):
    return np.array([50.0,100.0,-25.0,150.0,200.0,-50.0] + list([0.25]))

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

    self.elastic = elasticity.IsotropicLinearElasticModel(self.mu,
        "shear", self.K, "bulk")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model = models.SmallStrainRateIndependentPlasticity(self.elastic, flow)

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

    self.elastic = elasticity.IsotropicLinearElasticModel(self.mu,
        "shear", self.K, "bulk")

    surface = surfaces.IsoJ2()
    hrule = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model = models.SmallStrainRateIndependentPlasticity(self.elastic, flow)

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

    self.elastic = elasticity.IsotropicLinearElasticModel(self.mu,
        "shear", self.K, "bulk")

    surface = surfaces.IsoJ2()
    hrule = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model = models.SmallStrainRateIndependentPlasticity(self.elastic, flow)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

  def gen_hist(self):
    return np.array(range(1,8)) / 8.0

  def gen_x(self):
    return np.array(range(1,9)) / 9.0


class TestRIChabocheLinear(unittest.TestCase, CommonMatModel, CommonJacobian):
  """
    Test Chaboche with linear isotropic hardening
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
    self.cs = [10.0, 2.0]
    self.rs = [5.0, 1.0]
    self.gmodels = [hardening.ConstantGamma(g) for g in self.rs]
    self.As = [0.0] * self.n
    self.ns = [1.0] * self.n

    self.elastic = elasticity.IsotropicLinearElasticModel(self.mu,
        "shear", self.K, "bulk")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    hmodel = hardening.Chaboche(iso, self.cs, self.gmodels, 
        self.As, self.ns)

    flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hmodel)

    self.model = models.SmallStrainRateIndependentPlasticity(self.elastic, flow)

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

# Something funny with the Jacobian here
class TestCreepPlasticityJ2LinearPowerLaw(unittest.TestCase, CommonMatModel):
  """
    Test the combined creep/plasticity algorithm with J2 plasticity with
    isotropic hardening and power law creep 
  """
  def setUp(self):
    self.hist0 = np.zeros((13,))

    self.A = 1.85e-10
    self.n = 2.5

    self.smodel = creep.PowerLawCreep(self.A, self.n)
    self.cmodel = creep.J2CreepModel(self.smodel)

    self.E = 150000.0
    self.nu = 0.3
    self.sY = 200.0
    self.H = self.E / 50.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E, "youngs", 
        self.nu, "poissons")

    self.surface = surfaces.IsoJ2()
    self.iso = hardening.LinearIsotropicHardeningRule(self.sY, self.H)
    self.flow = ri_flow.RateIndependentAssociativeFlow(self.surface, self.iso)

    self.pmodel = models.SmallStrainRateIndependentPlasticity(self.elastic, 
        self.flow)

    self.model = models.SmallStrainCreepPlasticity(self.elastic, self.pmodel,
        self.cmodel)

    self.efinal = np.array([0.1,-0.05,0.02,-0.03,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

  def gen_hist(self):
    return np.array(range(1,7) + [1.0] + range(1,7)) / 7.0

  def gen_x(self):
    return np.array(range(1,7)) / 7.0

  def gen_start_stress(self):
    return np.zeros((6,)) + 100.0

  def gen_start_strain(self):
    return np.zeros((6,)) + 0.01

class TestCreepPlasticityPerfect(unittest.TestCase, CommonMatModel):
  """
    Test the combined creep/plasticity algorithm with J2 plasticity with
    perfect plasticity
  """
  def setUp(self):
    self.hist0 = np.zeros((6,))

    self.A = 1.85e-10
    self.n = 2.5

    self.smodel = creep.PowerLawCreep(self.A, self.n)
    self.cmodel = creep.J2CreepModel(self.smodel)

    self.E = 150000.0
    self.nu = 0.3
    self.sY = 200.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E, "youngs", 
        self.nu, "poissons")
    self.surface = surfaces.IsoJ2()

    self.pmodel = models.SmallStrainPerfectPlasticity(self.elastic, 
        self.surface, self.sY)

    self.model = models.SmallStrainCreepPlasticity(self.elastic, 
        self.pmodel, self.cmodel)

    self.efinal = np.array([0.1,-0.05,0.02,-0.05,0.1,-0.15])
    self.tfinal = 10.0
    self.T = 300.0
    self.nsteps = 10

  def gen_hist(self):
    return np.array(range(1,7) + [1.0] + range(1,7)) / 7.0

  def gen_x(self):
    return np.array(range(1,7)) / 7.0

  def gen_start_stress(self):
    return np.zeros((6,)) + 100.0

  def gen_start_strain(self):
    return np.zeros((6,)) + 0.01

class TestDirectIntegrateChaboche(unittest.TestCase, CommonMatModel, CommonJacobian):
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

    self.model = models.GeneralIntegrator(self.elastic, flow)

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

    self.E = 92000.0
    self.nu = 0.3

    self.mu = self.E/(2*(1+self.nu))
    self.K = self.E/(3*(1-2*self.nu))

    self.s0 = 180.0
    self.R = 150.0
    self.d = 10.0

    self.n = 2.0
    self.eta = 200.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.mu,
        "shear", self.K, "bulk")

    surface = surfaces.IsoJ2()
    hrule = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)
    g = visco_flow.GPowerLaw(self.n, self.eta)

    vmodel = visco_flow.PerzynaFlowRule(surface, hrule, g)

    flow = general_flow.TVPFlowRule(self.elastic, vmodel)
    self.model = models.GeneralIntegrator(self.elastic, flow)

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

    shear = interpolate.PolynomialInterpolate(mu_poly)
    bulk = interpolate.PolynomialInterpolate(K_poly)
    self.elastic = elasticity.IsotropicLinearElasticModel(shear, "shear", 
        bulk, "bulk")

    flow = general_flow.TVPFlowRule(self.elastic, vmodel)

    self.model = models.GeneralIntegrator(self.elastic, flow)

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

class TestKMSwitch(unittest.TestCase, CommonMatModel):
  """
    Test the model with a switch
  """
  def setUp(self):
    # Fully-defined perfectly plastic model
    Epoly = [-78.2759, 236951.0]
    nu = 0.3
    A = -9.6187
    B = -1.4819
    C = -5.0486
    g0 = 0.3708

    b = 0.248 * 1.0e-6
    kboltz = 1.38064e-23 * 1000.0
    eps0 = 1.0e10

    # Temperature range over which to consider (K)
    Tmin = 550.0
    Tmax = 950.0
    Trange = np.linspace(Tmin, Tmax)

    # Elastic
    E_m = interpolate.PolynomialInterpolate(Epoly)
    nu_m = interpolate.ConstantInterpolate(nu)
    elastic_m = elasticity.IsotropicLinearElasticModel(E_m, "youngs",
        nu_m, "poissons")
    self.elastic = elastic_m

    # Rate sensitivity interpolates values
    mu_values = np.array([elastic_m.G(T) for T in Trange])
    n_values = -mu_values*b**3.0 / (kboltz * Trange * A)
    eta_values = np.exp(B) * eps0 ** (kboltz * Trange * A / (mu_values * b**3.0)) * mu_values

    # Rate independent interpolate values
    flow_stress = mu_values * np.exp(C)
    
    # Common objects
    surface = surfaces.IsoKinJ2()
    hmodulus = interpolate.PolynomialInterpolate([-10.0, 12000.0])

    # Setup visco model
    n_interp = interpolate.PiecewiseLinearInterpolate(list(Trange), list(n_values))
    eta_interp = interpolate.PiecewiseLinearInterpolate(list(Trange), list(eta_values))
    eta_m = visco_flow.ConstantFluidity(eta_interp)

    iso_rd = hardening.LinearIsotropicHardeningRule(
        interpolate.ConstantInterpolate(0.0),
        hmodulus)
    hard_rd = hardening.Chaboche(iso_rd,
        [interpolate.ConstantInterpolate(0.0)], 
        [hardening.ConstantGamma(interpolate.ConstantInterpolate(0.0))],
        [interpolate.ConstantInterpolate(0.0)],
        [interpolate.ConstantInterpolate(1.0)])

    visco_rd = visco_flow.ChabocheFlowRule(surface, hard_rd, eta_m, n_interp) 
    general_rd = general_flow.TVPFlowRule(elastic_m, visco_rd)

    rate_dependent = models.GeneralIntegrator(elastic_m, general_rd)

    # Setup rate independent
    sy_interp = interpolate.PiecewiseLinearInterpolate(list(Trange), list(flow_stress))
    iso_ri = hardening.LinearIsotropicHardeningRule(sy_interp, hmodulus)
    hard_ri = hardening.Chaboche(iso_ri,
        [interpolate.ConstantInterpolate(0.0)], 
        [hardening.ConstantGamma(interpolate.ConstantInterpolate(0.0))],
        [interpolate.ConstantInterpolate(0.0)],
        [interpolate.ConstantInterpolate(1.0)])
    flow_ri = ri_flow.RateIndependentNonAssociativeHardening(surface, hard_ri)
    
    rate_independent = models.SmallStrainRateIndependentPlasticity(elastic_m,
        flow_ri)

    # Combined model
    self.model = models.KMRegimeModel(elastic_m, [rate_independent, rate_dependent],
        [g0], kboltz, b, eps0)

    self.efinal = np.array([0.05,0,0,0.02,0,-0.01])
    self.tfinal = 10.0
    self.T = 700.0
    self.nsteps = 100

  def gen_hist(self):
    h = np.array([40.0,20,-30,40.0,5.0,2.0,40.0])
    h[:6] = make_dev(h[:6])
    return h
