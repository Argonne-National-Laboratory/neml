import sys
sys.path.append('..')

from neml import elasticity, interpolate
from neml.math import tensors
import unittest

from common import *

import numpy as np
import numpy.linalg as la

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

  def test_tensor_C(self):
    CT = self.model.C_tensor(self.T)
    C = self.model.C(self.T)
    self.assertEqual(CT, tensors.SymSym(C))

  def test_tensor_S(self):
    ST = self.model.S_tensor(self.T)
    S = self.model.S(self.T)
    self.assertEqual(ST, tensors.SymSym(S))

class TestIsotropicConstantModel(CommonElasticity, unittest.TestCase):
  def setUp(self):
    self.mu = 29000.0
    self.K = 64000.0
    self.T = 325.0

    self.model = elasticity.IsotropicLinearElasticModel(self.mu, 
        "shear", self.K, "bulk")

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
        self.E, "youngs", self.nu, "poissons")
    self.model_GK = elasticity.IsotropicLinearElasticModel(
        self.mu, "shear", self.K, "bulk")

    self.T = 300.0

  def test_equivalent_C(self):
    C1 = self.model_Ev.C(self.T)
    C2 = self.model_GK.C(self.T)
    self.assertTrue(np.allclose(C1,C2))

  def test_equivalent_S(self):
    S1 = self.model_Ev.S(self.T)
    S2 = self.model_GK.S(self.T)
    self.assertTrue(np.allclose(S1,S2))
  
  def test_equivalent_modulii(self):
    self.assertTrue(np.isclose(self.E, self.model_Ev.E(self.T)))
    self.assertTrue(np.isclose(self.E, self.model_GK.E(self.T)))
    self.assertTrue(np.isclose(self.nu, self.model_Ev.nu(self.T)))
    self.assertTrue(np.isclose(self.nu, self.model_GK.nu(self.T)))
    self.assertTrue(np.isclose(self.mu, self.model_Ev.G(self.T)))
    self.assertTrue(np.isclose(self.mu, self.model_GK.G(self.T)))
    self.assertTrue(np.isclose(self.K, self.model_Ev.K(self.T)))
    self.assertTrue(np.isclose(self.K, self.model_GK.K(self.T)))
