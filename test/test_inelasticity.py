#!/usr/bin/env python3

from neml import history, interpolate
from neml.math import tensors, rotations
from neml.cp import crystallography, slipharden, sliprules, inelasticity

from common import differentiate
from nicediff import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonInelastic(object):
  def test_d_p_d_stress(self):
    nd = diff_symmetric_symmetric(lambda s: self.model.d_p(s, self.Q, self.H, self.L, self.T),
        self.S)
    d = self.model.d_d_p_d_stress(self.S, self.Q, self.H, self.L, self.T)
    self.assertEqual(nd, d)

  def test_d_p_d_history(self):
    nd = diff_symmetric_history(lambda h: self.model.d_p(self.S, self.Q, h, self.L, self.T),
        self.H)
    d = np.array(self.model.d_d_p_d_history(self.S, self.Q, self.H, self.L, self.T))
    
    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

  def test_d_hist_rate_d_stress(self):
    nd = diff_history_symmetric(lambda s: self.model.history_rate(s, self.Q, self.H, self.L, self.T),
        self.S)
    d = np.array(self.model.d_history_rate_d_stress(self.S, self.Q, self.H, self.L, self.T))

    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

  def test_d_hist_rate_d_hist(self):
    nd = diff_history_history(lambda h: self.model.history_rate(self.S, self.Q, h, self.L, self.T),
        self.H)
    d = np.array(self.model.d_history_rate_d_history(self.S, self.Q, self.H, self.L, self.T))

    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

class TestNoInelasticity(unittest.TestCase, CommonInelastic):
  def setUp(self):
    self.model = inelasticity.NoInelasticity()

    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))

    self.T = 300.0

    self.H = history.History()

  def test_d_p(self):
    self.assertEqual(tensors.Symmetric(np.zeros((3,3))),
        self.model.d_p(self.S, self.Q, self.H, self.L, self.T))

  def test_w_p(self):
    self.assertEqual(tensors.Skew(np.zeros((3,3))),
        self.model.w_p(self.S, self.Q, self.H, self.L, self.T))

  def test_hist_rate(self):
    h1 = history.History()
    h2 = self.model.history_rate(self.S, self.Q, self.H, self.L, self.T)

    self.assertTrue(np.allclose(np.array(h1), np.array(h2)))

class TestAsaroInelasticity(unittest.TestCase, CommonInelastic):
  def setUp(self):
    self.strength = 35.0
    self.H = history.History()
    self.H.add_scalar("strength")
    self.H.set_scalar("strength", self.strength)

    self.tau0 = 10.0
    self.tau_sat = 50.0
    self.b = 2.5

    self.strengthmodel = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0)
    
    self.g0 = 1.0
    self.n = 3.0
    self.slipmodel = sliprules.PowerLawSlipRule(self.strengthmodel, self.g0, self.n)

    self.model = inelasticity.AsaroInelasticity(self.slipmodel)

    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))

    self.T = 300.0

  def test_d_p(self):
    d = tensors.Symmetric(np.zeros((3,3)))
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        d += self.slipmodel.slip(g, i, self.S, self.Q, self.H, self.L, self.T) * self.L.M(g, i, self.Q)
    
    self.assertEqual(d, 
        self.model.d_p(self.S, self.Q, self.H, self.L, self.T))

  def test_w_p(self):
    w = tensors.Skew(np.zeros((3,3)))
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        w += self.slipmodel.slip(g, i, self.S, self.Q, self.H, self.L, self.T) * self.L.N(g, i, self.Q)
    
    self.assertEqual(w, 
        self.model.w_p(self.S, self.Q, self.H, self.L, self.T))

  def test_hist_rate(self):
    h1 = self.slipmodel.hist_rate(self.S, self.Q, self.H, self.L, self.T)
    h2 = self.model.history_rate(self.S, self.Q, self.H, self.L, self.T)

    self.assertTrue(np.allclose(np.array(h1), np.array(h2)))
