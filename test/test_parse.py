from neml import neml, elasticity, ri_flow, hardening, surfaces, parse

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

    for i,m in enumerate(np.linspace(0,1,self.nsteps)):
      t_np1 = self.tmax * m
      strain_np1 = self.emax * m
      
      stress_np11, hist_np11, A_np11 = self.model1.update_sd(strain_np1,
          strain_n, self.T, self.T, t_np1, t_n, stress_n1, hist_n1)

      stress_np12, hist_np12, A_np12 = self.model2.update_sd(strain_np1,
          strain_n, self.T, self.T, t_np1, t_n, stress_n2, hist_n2)

      self.assertTrue(np.allclose(stress_np11, stress_np12))
      self.assertTrue(np.allclose(hist_np11, hist_np12))
      self.assertTrue(np.allclose(A_np11, A_np12))
      
      stress_n1 = stress_np11
      hist_n1 = hist_np11
      stress_n2 = stress_np12
      hist_n2 = hist_np12

      strain_n = strain_np1


class TestJ2Iso(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_j2iso")
  
    mu = 40000.0
    K = 84000.0

    ys = 100.0
    H = 1000.0

    shear = elasticity.ConstantShearModulus(mu)
    bulk = elasticity.ConstantBulkModulus(K)
    elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

    surface = surfaces.IsoJ2()
    hrule = hardening.LinearIsotropicHardeningRule(ys, H)
    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.model2 = neml.SmallStrainRateIndependentPlasticity(elastic, flow)

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

    shear = elasticity.ConstantShearModulus(mu)
    bulk = elasticity.ConstantBulkModulus(K)
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

class TestJ2Combined(CompareMats, unittest.TestCase):
  def setUp(self):
    self.model1 = parse.parse_xml("test/examples.xml", "test_nonassri")
  
    mu = 40000.0
    K = 84000.0

    ys = 100.0
    r = 100.0
    d = 1000.0
    
    cs = [5.0,10.0]

    gs = [1000.0,1000.0]

    shear = elasticity.ConstantShearModulus(mu)
    bulk = elasticity.ConstantBulkModulus(K)
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


