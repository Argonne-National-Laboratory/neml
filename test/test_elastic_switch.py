import sys
sys.path.append("..")

from neml import models, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep, damage
from common import *

import unittest

class TestSimpleMaterials(unittest.TestCase):
  def setUp(self):
    self.G1 = 36000.0
    self.G2 = 25000.0

    self.K1 = 10000.0
    self.K2 = 12000.0

    self.elastic1 = elasticity.IsotropicLinearElasticModel(
        self.G1, "shear",
        self.K1, "bulk")
    self.elastic2 = elasticity.IsotropicLinearElasticModel(
        self.G2, "shear",
        self.K2, "bulk")
    
  def are_equal(self, e1, e2, T = 30):
    self.assertTrue(np.allclose(e1.C(T), e2.C(T)))

  def test_linear_elastic(self):
    model = models.SmallStrainElasticity(self.elastic1)
    self.are_equal(self.elastic1, model.elastic)
    model.set_elastic_model(self.elastic2)
    self.are_equal(self.elastic2, model.elastic)

  def test_perfect_plasicity(self):
    surface = surfaces.IsoJ2()
    sy = 100.0
    model = models.SmallStrainPerfectPlasticity(self.elastic1, surface, sy)
    self.are_equal(self.elastic1, model.elastic)
    model.set_elastic_model(self.elastic2)
    self.are_equal(self.elastic2, model.elastic)

  def test_ri_plasticity(self):
    surface = surfaces.IsoJ2()
    sy = 100.0
    K = 1000.0
    H = 500.0
    hrule = hardening.CombinedHardeningRule(
        hardening.LinearIsotropicHardeningRule(sy, K),
        hardening.LinearKinematicHardeningRule(H))
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)
    model = models.SmallStrainRateIndependentPlasticity(self.elastic1, 
        flow)
    self.are_equal(self.elastic1, model.elastic)
    model.set_elastic_model(self.elastic2)
    self.are_equal(self.elastic2, model.elastic)

class TestCompoundMaterials(unittest.TestCase):
  def setUp(self):
    self.G1 = 36000.0
    self.G2 = 25000.0

    self.K1 = 10000.0
    self.K2 = 12000.0

    self.elastic1 = elasticity.IsotropicLinearElasticModel(
        self.G1, "shear",
        self.K1, "bulk")
    self.elastic2 = elasticity.IsotropicLinearElasticModel(
        self.G2, "shear",
        self.K2, "bulk")

    self.emodel1 = models.SmallStrainElasticity(self.elastic1)
    self.emodel2 = models.SmallStrainElasticity(self.elastic2)

    self.tiny_strain = np.array([1e-6,0,0,0,0,0])
    self.T = 300.0
    self.t = 1e-9

  def very_close(self, model1, model2):
    self.assertTrue(np.allclose(self.small_step(model1), 
      self.small_step(model2)))

  def small_step(self, model):
    h0 = model.init_store()
    s0 = np.zeros((6,))
    e0 = np.zeros((6,))
    t0 = 0.0
    u0 = 0.0
    p0 = 0.0

    stress, hist, A, u, p = model.update_sd(
        self.tiny_strain, e0, self.T, self.T, self.t, t0,
        s0, h0, u0, p0)

    return stress
  
  def test_comparison(self):
    self.very_close(self.emodel1, self.emodel1)

  def test_creep_plasticity(self):
    surface = surfaces.IsoJ2()
    sy = 100.0
    bmodel = models.SmallStrainPerfectPlasticity(self.elastic1, surface, sy)
    A = 1.85e-10
    n = 12.0
    
    smodel = creep.PowerLawCreep(A, n)
    cmodel = creep.J2CreepModel(smodel)
    model = models.SmallStrainCreepPlasticity(self.elastic1, bmodel, cmodel)

    self.very_close(model, self.emodel1)

    model.set_elastic_model(self.elastic2)

    self.very_close(model, self.emodel2)

  def test_rd(self):
    s0 = 10.0
    R = 150.0
    d = 10.0

    n = 3.0
    eta = 200.0

    surface = surfaces.IsoJ2()
    hrule = hardening.VoceIsotropicHardeningRule(s0, R, d)
    g = visco_flow.GPowerLaw(n, eta)

    vmodel = visco_flow.PerzynaFlowRule(surface, hrule, g)
    flow = general_flow.TVPFlowRule(self.elastic1, vmodel)
    model = models.GeneralIntegrator(self.elastic1, flow)

    self.very_close(model, self.emodel1)
    
    model.set_elastic_model(self.elastic2)

    self.very_close(model, self.emodel2)

  def test_damage(self):
    s0 = 180.0
    Kp = 1000.0
    H = 1000.0

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
    kin = hardening.LinearKinematicHardeningRule(H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    bmodel = models.SmallStrainRateIndependentPlasticity(self.elastic1, 
        flow)

    W0 = 10.0
    k0 = 0.0001
    a0 = 2.0

    model1 = damage.NEMLExponentialWorkDamagedModel_sd(
        self.elastic1, W0, k0, 
        a0, bmodel)

    W02 = 10.0
    k02 = 0.001
    a02 = 1.5

    model2 = damage.NEMLExponentialWorkDamagedModel_sd(
        self.elastic1, W02, k02, 
        a02, bmodel)

    model = damage.CombinedDamageModel_sd(self.elastic1, 
        [model1, model2], bmodel)
    
    self.very_close(model, self.emodel1)
    
    model.set_elastic_model(self.elastic2)

    self.very_close(model, self.emodel2)
