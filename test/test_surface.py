import sys
sys.path.append('..')

from neml.surface import *
from common import *

import unittest
import numpy as np

class TestLinearKinIsoJ2(unittest.TestCase):
  def setUp(self):
    self.K0 = 450.0
    self.Kb = 10000.0
    self.Hb = 1000.0

    self.model = LinearKinIsoJ2(self.K0, self.Kb, self.Hb)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.K0, self.K0))
    self.assertTrue(np.isclose(self.model.Kb, self.Kb))
    self.assertTrue(np.isclose(self.model.Hb, self.Hb))
