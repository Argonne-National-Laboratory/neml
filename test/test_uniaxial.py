import sys
sys.path.append('..')

from neml import interpolate, solvers, models, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep, uniaxial
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonUniaxial(object):
  """
    Common uniaxial model tests
  """
  def test_tangent(self):
    hn = self.umodel.init_store()
    en = 0.0
    sn = 0.0
    un = 0.0
    pn = 0.0
    tn = 0.0
    Tn = self.T0

    n = 100
    emax = 0.1

    for i in range(self.nsteps):
      enp1 = en + self.de
      Tnp1 = Tn + self.dT
      tnp1 = tn + self.dt
      snp1, hnp1, Anp1, unp1, pnp1 = self.umodel.update(enp1, en, Tnp1, Tn, 
          tnp1, tn, sn, hn, un, pn)
      
      sfn = lambda e: self.umodel.update(e, en, Tnp1, Tn, tnp1, tn,
          sn, hn, un, pn)[0]
      nA = differentiate(sfn, enp1)
      
      self.assertTrue(np.isclose(nA, Anp1, rtol = 1.0e-3))

      sn = snp1
      hn = np.copy(hnp1)
      en = enp1
      Tn = Tnp1
      tn = tnp1
      un = unp1
      pn = pnp1

class TestUniaxialRI(CommonUniaxial, unittest.TestCase):
  """
    Test with a simple uniaxial model
  """
  def setUp(self):
    E = 200000.0
    nu = 0.27

    mu = E / (2 * (1.0 + nu))
    K = E / (3 * (1 - 2 * nu))

    s0 = 300.0
    Kp = 0.0
    c = [30000.0]
    r = [60.0]
    A = [0.0]
    n = [1.0]
    
    elastic = elasticity.IsotropicLinearElasticModel(mu, "shear",
        K, "bulk")
    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
    gmodels = [hardening.ConstantGamma(g) for g in r]
    hrule = hardening.Chaboche(iso, c, gmodels, A, n)

    flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hrule)
    self.model = models.SmallStrainRateIndependentPlasticity(elastic, flow, verbose = False)
    self.umodel = uniaxial.UniaxialModel(self.model)
    
    self.de = 0.001
    self.dt = 1.0
    self.dT = 0.0
    self.T0 = 0.0
    self.nsteps = 100
