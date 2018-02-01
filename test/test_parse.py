from neml import solvers, interpolate, neml, elasticity, ri_flow, hardening, surfaces, parse, visco_flow, general_flow, creep, damage

import unittest
import numpy as np

class CompareMats(object):
  def test_same(self):
    t_n = 0.0
    stress_n1 = np.zeros((6,))
    hist_n1 = self.model1.init_store()

    stress_n2 = np.zeros((6,))
    hist_n2 = self.model2.init_store()

    strain_n = np.zeros((6,))

    u_n1 = 0.0
    u_n2 = 0.0
    p_n1 = 0.0
    p_n2 = 0.0

    for i,m in enumerate(np.linspace(0,1,self.nsteps)):
      t_np1 = self.tmax * m
      strain_np1 = self.emax * m
      
      stress_np11, hist_np11, A_np11, u_np11, p_np11 = self.model1.update_sd(
          strain_np1, strain_n, self.T, self.T, t_np1, t_n, stress_n1, hist_n1,
          u_n1, p_n1)

      stress_np12, hist_np12, A_np12, u_np12, p_np12 = self.model2.update_sd(
          strain_np1, strain_n, self.T, self.T, t_np1, t_n, stress_n2, hist_n2,
          u_n2, p_n2)

      self.assertTrue(np.allclose(stress_np11, stress_np12))
      self.assertTrue(np.allclose(hist_np11, hist_np12))
      self.assertTrue(np.allclose(A_np11, A_np12))
      self.assertTrue(np.isclose(u_np11, u_np12))
      self.assertTrue(np.isclose(p_np11, p_np12))
      
      stress_n1 = stress_np11
      hist_n1 = hist_np11
      stress_n2 = stress_np12
      hist_n2 = hist_np12

      u_n1 = u_np11
      u_n2 = u_np12
      p_n1 = p_np11
      p_n2 = p_np12

      strain_n = strain_np1

class TestMarkDamage(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_markdamage")

    E = 92000.0
    nu = 0.3

    s0 = 180.0
    Kp = 1000.0

    E = elasticity.YoungsModulus(E)
    v = elasticity.PoissonsRatio(nu)
    elastic = elasticity.IsotropicLinearElasticModel(E, v)

    surface = surfaces.IsoJ2()
    hrule = hardening.LinearIsotropicHardeningRule(s0, Kp)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    bmodel = neml.SmallStrainRateIndependentPlasticity(elastic, 
        flow)

    C = 8.0e-6
    n = 2.2
    m = 1.0
    alpha = 3.0
    beta = 2.0
    r0 = 0.01
    self.model2 = damage.MarkFatigueDamageModel_sd(C, m, n, alpha, beta, r0, 
        bmodel)

    self.T = 300.0
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.05,0,0,0,0,0])


class TestPowerLawDamage(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_powerdamage")

    E = 92000.0
    nu = 0.3

    s0 = 180.0
    Kp = 1000.0

    E = elasticity.YoungsModulus(E)
    v = elasticity.PoissonsRatio(nu)
    elastic = elasticity.IsotropicLinearElasticModel(E, v)

    surface = surfaces.IsoJ2()
    hrule = hardening.LinearIsotropicHardeningRule(s0, Kp)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    bmodel = neml.SmallStrainRateIndependentPlasticity(elastic, 
        flow)

    a = 2.2
    A = 2e-5
    self.model2 = damage.NEMLPowerLawDamagedModel_sd(A, a, bmodel, elastic)

    self.T = 300.0
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.05,0,0,0,0,0])

class TestJ2Iso(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_j2iso")
  
    mu = 40000.0
    K = 84000.0

    ys = 100.0
    H = 1000.0

    shear = elasticity.ShearModulus(mu)
    bulk = elasticity.BulkModulus(K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoJ2()
    hrule = hardening.LinearIsotropicHardeningRule(ys, H)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model2 = neml.SmallStrainRateIndependentPlasticity(elastic, flow)

    self.T = 300.0
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.1,0,0,0,0,0])

class TestJ2Iso(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_j2isocomb")
  
    E = 150000.0
    nu = 0.3

    ys = 100.0
    H = 100.0
    r = 100.0
    d = 1000.0

    ym = elasticity.YoungsModulus(E)
    pr = elasticity.PoissonsRatio(nu)
    elastic = elasticity.IsotropicLinearElasticModel(ym, pr)

    surface = surfaces.IsoJ2()
    hrule1 = hardening.LinearIsotropicHardeningRule(ys, H)
    hrule2 = hardening.VoceIsotropicHardeningRule(0.0, r, d)
    hrule = hardening.CombinedIsotropicHardeningRule([hrule1,hrule2])
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model2 = neml.SmallStrainRateIndependentPlasticity(elastic, flow)

    self.T = 300.0
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.1,0,0,0,0,0])

class TestCreepPlasticity(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_creep_plasticity")
  
    A = 1.85e-10
    n = 2.5

    smodel = creep.PowerLawCreep(A, n)
    cmodel = creep.J2CreepModel(smodel)

    E = 150000.0
    nu = 0.3
    sY = 200.0
    H = E / 50.0

    youngs = elasticity.YoungsModulus(E)
    poisson = elasticity.PoissonsRatio(nu)
    elastic = elasticity.IsotropicLinearElasticModel(youngs, poisson)
    surface = surfaces.IsoJ2()
    iso = hardening.LinearIsotropicHardeningRule(sY, H)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, iso)

    pmodel = neml.SmallStrainRateIndependentPlasticity(elastic, flow)
    
    self.model2 = neml.SmallStrainCreepPlasticity(pmodel, cmodel)

    self.T = 300.0
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.1,0,0,0,0,0])

class TestJ2Combined(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_j2comb")
  
    mu = 40000.0
    K = 84000.0

    ys = 100.0
    r = 100.0
    d = 1000.0

    KH = 1000.0

    shear = elasticity.ShearModulus(mu)
    bulk = elasticity.BulkModulus(K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoKinJ2()
    iso = hardening.VoceIsotropicHardeningRule(ys, r, d)
    kin = hardening.LinearKinematicHardeningRule(KH)
    hrule = hardening.CombinedHardeningRule(iso, kin)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model2 = neml.SmallStrainRateIndependentPlasticity(elastic, flow)

    self.T = 300.0
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.1,0,0,0,0,0])

class TestRIChaboche(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_nonassri")
  
    mu = 40000.0
    K = 84000.0

    ys = 100.0
    r = 100.0
    d = 1000.0
    
    cs = [5.0,10.0]

    gs = [1000.0,1000.0]

    shear = elasticity.ShearModulus(mu)
    bulk = elasticity.BulkModulus(K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoKinJ2()
    iso = hardening.VoceIsotropicHardeningRule(ys, r, d)
    
    gammas = [hardening.ConstantGamma(g) for g in gs]
    hmodel = hardening.Chaboche(iso, cs, gammas)
    flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hmodel)

    self.model2 = neml.SmallStrainRateIndependentPlasticity(elastic, flow)

    self.T = 300.0
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.1,0,0,0,0,0])

class TestYaguchi(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_yaguchi")
  
    mu_poly = [-0.11834615, 115.5, 48807.69]
    K_poly = [-0.256417, 250.25, 105750.0]

    shear = elasticity.ShearModulus(interpolate.PolynomialInterpolate(mu_poly))
    bulk = elasticity.BulkModulus(interpolate.PolynomialInterpolate(K_poly))
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)
    
    vmodel = visco_flow.YaguchiGr91FlowRule()

    flow = general_flow.TVPFlowRule(elastic, vmodel)

    self.model2 = neml.GeneralIntegrator(flow)

    self.T = 500.0
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.1,0,0,0,0,0])

class TestRDChaboche(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_rd_chaboche")
  
    mu = 60384.61
    K = 130833.3

    shear = elasticity.ShearModulus(mu)
    bulk = elasticity.BulkModulus(K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    r = -80.0
    d = 3.0

    Cs = [135.0e3, 61.0e3, 11.0e3]
    gs = [5.0e4, 1100.0, 1.0]
    gmodels = [hardening.ConstantGamma(g) for g in gs]

    eta = 701.0

    n = 10.5
    
    surface = surfaces.IsoKinJ2()
    iso = hardening.VoceIsotropicHardeningRule(0.0, r, d)

    hmodel = hardening.Chaboche(iso, Cs, gmodels)

    fluidity = visco_flow.ConstantFluidity(eta)
    
    vmodel = visco_flow.ChabocheFlowRule(surface, hmodel, fluidity, n)
    flow = general_flow.TVPFlowRule(elastic, vmodel)

    self.model2 = neml.GeneralIntegrator(flow)

    self.T = 550.0 + 273.15
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.1,0,0,0,0,0])


class TestPerzyna(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_perzyna")
  
    mu = 40000.0
    K = 84000.0

    shear = elasticity.ShearModulus(mu)
    bulk = elasticity.BulkModulus(K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)
    
    sy = 100.0
    r = 100.0
    d = 1000.0
    KK = 1000.0

    n = 5.0
    eta = 500.0
    
    surface = surfaces.IsoKinJ2()
    iso = hardening.VoceIsotropicHardeningRule(sy, r, d)
    kin = hardening.LinearKinematicHardeningRule(KK)

    hmodel = hardening.CombinedHardeningRule(iso, kin)

    gmodel = visco_flow.GPowerLaw(n)
    
    vmodel = visco_flow.PerzynaFlowRule(surface, hmodel, gmodel, eta)
    flow = general_flow.TVPFlowRule(elastic, vmodel)

    self.model2 = neml.GeneralIntegrator(flow)

    self.T = 550.0 + 273.15
    self.tmax = 10.0
    self.nsteps = 100.0
    self.emax = np.array([0.1,0,0,0,0,0])


class TestPerfect(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_perfect")

    E = [-100, 100000]
    nu = 0.3

    youngs = elasticity.YoungsModulus(interpolate.PolynomialInterpolate(E))
    poissons = elasticity.PoissonsRatio(nu)
    elastic = elasticity.IsotropicLinearElasticModel(youngs, poissons)

    surface = surfaces.IsoJ2()

    Ts = [100.0, 300.0, 500.0, 700.0]
    Sys = [1000.0, 120.0, 60.0, 30.0]

    yields = interpolate.PiecewiseLinearInterpolate(Ts, Sys)

    self.model2 = neml.SmallStrainPerfectPlasticity(elastic, surface, yields)

    self.T = 550.0
    self.tmax = 10.0
    self.nsteps = 50.0
    self.emax = np.array([0.1,0.05,0,-0.025,0,0])

  def test_alpha(self):
    self.assertTrue(np.isclose(self.model1.alpha(self.T), 0.1))

class TestPerfectCreep(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_pcreep")

    E = [-100, 100000]
    nu = 0.3

    youngs = elasticity.YoungsModulus(interpolate.PolynomialInterpolate(E))
    poissons = elasticity.PoissonsRatio(nu)
    elastic = elasticity.IsotropicLinearElasticModel(youngs, poissons)

    surface = surfaces.IsoJ2()

    Ts = [100.0, 300.0, 500.0, 700.0]
    Sys = [1000.0, 120.0, 60.0, 30.0]

    yields = interpolate.PiecewiseLinearInterpolate(Ts, Sys)

    pmodel = neml.SmallStrainPerfectPlasticity(elastic, surface, yields)

    self.T = 550.0
    self.tmax = 10.0
    self.nsteps = 50.0
    self.emax = np.array([0.1,0.05,0,-0.025,0,0])

    A = 1.85e-10
    n = 2.5

    smodel = creep.PowerLawCreep(A, n)
    cmodel = creep.J2CreepModel(smodel)

    self.model2 = neml.SmallStrainCreepPlasticity(pmodel, cmodel)
