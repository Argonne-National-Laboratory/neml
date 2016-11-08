import sys
sys.path.append('..')

from neml import elasticity
import unittest

from common import *

import numpy as np
import numpy.linalg as la

class TestConstantShearModulus(unittest.TestCase):
  def setUp(self):
    self.mu = 29000.0
    self.T = 325.0
    self.model = elasticity.ConstantShearModulus(self.mu)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.mu, self.mu))

  def test_modulus(self):
    self.assertTrue(np.isclose(self.model.mu, self.model.modulus(self.T)))

class TestConstantBulkModulus(unittest.TestCase):
  def setUp(self):
    self.K = 64000.0
    self.T = 325.0
    self.model = elasticity.ConstantBulkModulus(self.K)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.K, self.K))

  def test_modulus(self):
    self.assertTrue(np.isclose(self.model.K, self.model.modulus(self.T)))

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

    self.shear = elasticity.ConstantShearModulus(self.mu)
    self.bulk = elasticity.ConstantBulkModulus(self.K)

    self.model = elasticity.IsotropicLinearElasticModel(self.shear, self.bulk)

  def test_modulii(self):
    S = self.model.S(self.T)
    E = 9*self.K*self.mu/(3*self.K + self.mu)
    nu = (3*self.K - 2*self.mu)/(2*(3*self.K+self.mu))

    self.assertTrue(np.isclose(S[0,0], 1/E))
    self.assertTrue(np.isclose(S[0,1], -nu/E))
    self.assertTrue(np.isclose(S[3,3], (1+nu)/E))
