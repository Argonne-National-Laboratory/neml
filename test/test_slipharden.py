#!/usr/bin/env python3

from neml import history, interpolate
from neml.math import tensors, rotations
from neml.cp import crystallography, slipharden, sliprules

from common import differentiate
from nicediff import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonPlasticSlipHardening():
  def sum_slip(self):
    return sum(self.sliprule.slip(g, i, self.S, self.Q, self.H, self.L, 
      self.T) for g in range(self.L.ngroup) for i in range(self.L.nslip(g)))

  def test_hist_rate(self):
    ss = self.sum_slip()

    self.assertTrue(np.isclose(self.model.hist_rate(self.S, self.Q, self.H,
      self.L, self.T, self.sliprule), self.model.hist_factor(self.H.get_scalar("strength"), 
        self.L, self.T) * ss))

  def test_d_hist_rate_d_stress(self):
    dfn = lambda s: self.model.hist_rate(s, self.Q, self.H,
        self.L, self.T, self.sliprule)
    nd = diff_scalar_symmetric(dfn, self.S)
    d = self.model.d_hist_rate_d_stress(self.S, self.Q, self.H, self.L,
        self.T, self.sliprule)
    self.assertEqual(d,nd)

  def test_d_hist_rate_d_hist(self):
    dfn = lambda h: self.model.hist_rate(self.S, self.Q, h, 
        self.L, self.T, self.sliprule)
    nd = diff_scalar_history(dfn, self.H)
    
    d = self.model.d_hist_rate_d_strength(self.S, self.Q, self.H, self.L,
        self.T, self.sliprule)
    
    print(np.array(d))
    print(np.array(nd))
    self.assertTrue(False)

class TestVoceHardening(unittest.TestCase, CommonPlasticSlipHardening):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
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

  def test_initialize_hist(self):
    H = history.History()
    self.model.populate_history(H)
    self.model.init_history(H)
    self.assertTrue(np.isclose(H.get_scalar('strength'), 0.0))
  
  def test_factor(self):
    self.assertTrue(np.isclose(
      self.model.hist_factor(self.strength, self.L, self.T),
      self.b * (self.tau_sat - self.strength) + self.tau0))

  def test_d_factor(self):
    nd = differentiate(lambda s: self.model.hist_factor(s, self.L, self.T), 
        self.strength)
    d = self.model.d_hist_factor(self.strength, self.L, self.T)

    self.assertTrue(np.isclose(nd, d))

