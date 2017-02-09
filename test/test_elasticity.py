import sys
sys.path.append('..')

from neml import elasticity, interpolate
import unittest

from common import *

import numpy as np
import numpy.linalg as la

class TestConstantShearModulus(unittest.TestCase):
  def setUp(self):
    self.mu = 29000.0
    self.T = 325.0
    self.model = elasticity.ShearModulus(self.mu)

  def test_modulus(self):
    self.assertTrue(np.isclose(self.mu, self.model.modulus(self.T)))

class TestPolyShearModulus(unittest.TestCase):
  def setUp(self):
    self.coefs = [3.077e-1, 3.003e2, 1.269e5]
    self.T = 525.0
    self.model = elasticity.ShearModulus(interpolate.PolynomialInterpolate(self.coefs))

  def test_modulus(self):
    self.assertTrue(np.isclose(self.model.modulus(self.T), 
      np.polyval(self.coefs, self.T)))

class TestConstantBulkModulus(unittest.TestCase):
  def setUp(self):
    self.K = 64000.0
    self.T = 325.0
    self.model = elasticity.BulkModulus(self.K)

  def test_modulus(self):
    self.assertTrue(np.isclose(self.K, self.model.modulus(self.T)))

class TestPolyBulkModulus(unittest.TestCase):
  def setUp(self):
    self.coefs = [3.077e-2, 3.003e3, 1.269e2]
    self.T = 525.0
    self.model = elasticity.BulkModulus(interpolate.PolynomialInterpolate(self.coefs))

  def test_modulus(self):
    self.assertTrue(np.isclose(self.model.modulus(self.T), 
      np.polyval(self.coefs, self.T)))

class CommonElasticity(object):
  """
    Tests that could apply to any elastic model.
  """
  def test_C2S(self):
    C = self.model.C(self.T)
    self.assertTrue(np.allclose(self.model.S(self.T), la.inv(C)))

  def test_S2C(self):
    S = self.model.S(self.T)
    self.assertTrue(np.allclose(self.model.C(self.T), la.inv(S)))

class TestIsotropicConstantModel(CommonElasticity, unittest.TestCase):
  def setUp(self):
    self.mu = 29000.0
    self.K = 64000.0
    self.T = 325.0

    self.shear = elasticity.ShearModulus(self.mu)
    self.bulk = elasticity.BulkModulus(self.K)

    self.model = elasticity.IsotropicLinearElasticModel(self.shear, self.bulk)

  def test_modulii(self):
    S = self.model.S(self.T)
    E = 9*self.K*self.mu/(3*self.K + self.mu)
    nu = (3*self.K - 2*self.mu)/(2*(3*self.K+self.mu))

    self.assertTrue(np.isclose(S[0,0], 1/E))
    self.assertTrue(np.isclose(S[0,1], -nu/E))
    self.assertTrue(np.isclose(S[3,3], (1+nu)/E))


class TestEquivalentDefinitions(unittest.TestCase):
  def setUp(self):
    self.E = 100000.0
    self.nu = 0.3

    self.mu = self.E / (2 * (1 + self.nu))
    self.K = self.E / (3 * (1 - 2.0 * self.nu))
    
    self.model_Ev = elasticity.IsotropicLinearElasticModel(
        elasticity.YoungsModulus(self.E), elasticity.PoissonsRatio(self.nu))
    self.model_GK = elasticity.IsotropicLinearElasticModel(
        elasticity.ShearModulus(self.mu), elasticity.BulkModulus(self.K))

    self.T = 300.0

  def test_equivalent_C(self):
    C1 = self.model_Ev.C(self.T)
    C2 = self.model_GK.C(self.T)
    self.assertTrue(np.allclose(C1,C2))

  def test_equivalent_S(self):
    S1 = self.model_Ev.S(self.T)
    S2 = self.model_GK.S(self.T)
    self.assertTrue(np.allclose(S1,S2))
