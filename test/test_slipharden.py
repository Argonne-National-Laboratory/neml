#!/usr/bin/env python3

from neml import history, interpolate
from neml.math import tensors, rotations
from neml.cp import crystallography, slipharden, sliprules

from common import differentiate

import unittest
import numpy as np
import numpy.linalg as la

class CommonPlasticSlipHardening():
  pass

class TestVoceHardening(unittest.TestCase):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))
    self.strength = 35.0
    self.H = history.History()
    self.H.add_scalar("strength")
    self.H.set_scalar("strength", self.strength)

    self.T = 300.0

    self.tau0 = 10.0
    self.tau_sat = 50.0
    self.b = 2.5

    self.model = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0)
    
    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)
  
  def test_factor(self):
    self.assertTrue(np.isclose(
      self.model.hist_factor(self.strength, self.L, self.T),
      self.b * (self.tau_sat - self.strength) + self.tau0))

  def test_d_factor(self):
    nd = differentiate(lambda s: self.model.hist_factor(s, self.L, self.T), 
        self.strength)
    d = self.model.d_hist_factor(self.strength, self.L, self.T)

    self.assertTrue(np.isclose(nd, d))

