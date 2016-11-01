import sys
sys.path.append('..')

from neml.shear import *
from common import *

import unittest
import numpy as np

class TestConstantShearModulus(unittest.TestCase):
  def setUp(self):
    self.mu = 15000.0
    self.T = 301.0

    self.model = ConstantShearModulus(self.mu)

  def test_properties(self):
    self.assertTrue(np.isclose(self.mu, self.model.mu))

  def test_modulus(self):
    self.assertTrue(np.isclose(self.mu, self.model.modulus(self.T)))
