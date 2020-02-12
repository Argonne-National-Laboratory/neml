#!/usr/bin/env python3

from neml import history, interpolate
from neml.math import tensors, rotations
from neml.cp import crystallography, slipharden, sliprules

from common import differentiate
from nicediff import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonSlipHardening():
  def test_d_hist_d_stress(self):
    d = np.array(self.model.d_hist_d_s(self.S, self.Q, self.H, self.L, self.T, self.sliprule))
    nd = diff_history_symmetric(lambda s: self.model.hist(s, self.Q, self.H, self.L, self.T,
      self.sliprule), self.S)
    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

  def test_d_hist_d_hist(self):
    d = np.array(self.model.d_hist_d_h(self.S, self.Q, self.H, self.L, self.T, self.sliprule))
    nd = diff_history_history(lambda h: self.model.hist(self.S, self.Q, h, self.L, self.T,
      self.sliprule), self.H)
    
    print(d)
    print(nd)
    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

  def test_d_hist_to_tau_d_hist(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        nd = diff_history_scalar(lambda h: self.model.hist_to_tau(g, i, h, self.T), self.H)
        d = self.model.d_hist_to_tau(g, i, self.H, self.T)
        self.assertTrue(np.allclose(np.array(nd), np.array(d)))

class CommonSlipSingleHardening():
  def test_d_hist_map(self):
    nd = diff_history_scalar(lambda h: self.model.hist_map(h, self.T), self.H)
    H = self.model.d_hist_map(self.H, self.T)
    self.assertTrue(np.allclose(np.array(H), np.array(nd)))

  def test_hist_to_tau(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        self.assertTrue(np.isclose(self.strength + self.static, 
          self.model.hist_to_tau(g, i, self.H, self.T)))

class CommonSlipSingleStrengthHardening():
  def test_hist_map(self):
    self.assertTrue(np.isclose(self.model.hist_map(self.H, self.T), self.strength + self.static))

  def test_hist(self):
    h = self.model.hist(self.S, self.Q, self.H, self.L, self.T, self.sliprule)
    self.assertTrue(np.isclose(h.get_scalar(self.vname), self.model.hist_rate(
      self.S, self.Q, self.H, self.L, self.T, self.sliprule)))

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
    nd = diff_history_scalar(dfn, self.H)
    
    d = self.model.d_hist_rate_d_hist(self.S, self.Q, self.H, self.L,
        self.T, self.sliprule)
    
    self.assertTrue(np.isclose(d.get_scalar(self.vname), nd.get_scalar(self.vname)))

class CommonPlasticSlipHardening():
  def sum_slip(self):
    return sum(np.abs(self.sliprule.slip(g, i, self.S, self.Q, self.H, self.L, 
      self.T)) for g in range(self.L.ngroup) for i in range(self.L.nslip(g)))

  def test_hist_rate(self):
    ss = self.sum_slip()

    self.assertTrue(np.isclose(self.model.hist_rate(self.S, self.Q, self.H,
      self.L, self.T, self.sliprule), self.model.hist_factor(self.H.get_scalar(self.vname), 
        self.L, self.T) * ss))

  def test_d_factor(self):
    nd = differentiate(lambda s: self.model.hist_factor(s, self.L, self.T), 
        self.strength)
    d = self.model.d_hist_factor(self.strength, self.L, self.T)

    self.assertTrue(np.isclose(nd, d))
    
class TestVoceHardening(unittest.TestCase, CommonPlasticSlipHardening, 
    CommonSlipSingleStrengthHardening, CommonSlipSingleHardening, CommonSlipHardening):
  def setUp(self):
    self.vname = "strength"

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

    self.static = self.tau0

    self.model = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0)
    
    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)

  def test_initialize_hist(self):
    H = history.History()
    self.model.populate_history(H)
    self.model.init_history(H)
    self.assertTrue(np.isclose(H.get_scalar('strength'), 0.0))

  def test_static_strength(self):
    self.assertTrue(np.isclose(self.model.static_strength(self.T), self.tau0))
  
  def test_factor(self):
    self.assertTrue(np.isclose(
      self.model.hist_factor(self.strength, self.L, self.T),
      self.b * (self.tau_sat - self.strength)))

class TestVoceHardeningMoveVar(unittest.TestCase, CommonPlasticSlipHardening, 
    CommonSlipSingleStrengthHardening, CommonSlipSingleHardening, CommonSlipHardening):
  def setUp(self):
    self.vname = "wee"
    
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))
    self.strength = 35.0
    self.H = history.History()
    self.H.add_scalar(self.vname)
    self.H.set_scalar(self.vname, self.strength)

    self.T = 300.0

    self.tau0 = 10.0
    self.tau_sat = 50.0
    self.b = 2.5

    self.static = self.tau0

    self.model = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0)
    
    self.model.set_varnames([self.vname])
    
    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)

  def test_get_vnames(self):
    self.assertEqual(self.model.varnames, [self.vname])

  def test_initialize_hist(self):
    H = history.History()
    self.model.populate_history(H)
    self.model.init_history(H)
    self.assertTrue(np.isclose(H.get_scalar(self.vname), 0.0))

  def test_static_strength(self):
    self.assertTrue(np.isclose(self.model.static_strength(self.T), self.tau0))
  
  def test_factor(self):
    self.assertTrue(np.isclose(
      self.model.hist_factor(self.strength, self.L, self.T),
      self.b * (self.tau_sat - self.strength)))

class TestMultiVoceHardening(unittest.TestCase, 
    CommonSlipSingleHardening, CommonSlipHardening):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))
    self.strength_0 = 35.0
    self.strength_1 = -5.0
    self.H = history.History()
    self.H.add_scalar("strength0")
    self.H.set_scalar("strength0", self.strength_0)
    self.H.add_scalar("strength1")
    self.H.set_scalar("strength1", self.strength_1)

    self.T = 300.0

    self.tau0_1 = 10.0
    self.tau_sat_1 = 50.0
    self.b_1 = 2.5

    self.tau0_2 = 5.0
    self.tau_sat_2 = -25.0
    self.b_2 = 1.5

    self.static = self.tau0_1 + self.tau0_2
    self.strength = self.strength_0 + self.strength_1

    self.model1 = slipharden.VoceSlipHardening(self.tau_sat_1, self.b_1, self.tau0_1)
    self.model2 = slipharden.VoceSlipHardening(self.tau_sat_2, self.b_2, self.tau0_2)

    self.model = slipharden.SumSlipSingleStrengthHardening([self.model1,
      self.model2])
    
    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)

  def test_initialize_hist(self):
    H = history.History()
    self.model.populate_history(H)
    self.model.init_history(H)
    self.assertTrue(np.isclose(H.get_scalar('strength0'), 0.0))
    self.assertTrue(np.isclose(H.get_scalar('strength1'), 0.0))

  def test_hist_def(self):
    hrate_1 = self.b_1 * (self.tau_sat_1 - self.strength_0) * self.sliprule.sum_slip(self.S, self.Q, self.H, self.L, self.T)
    hrate_2 = self.b_2 * (self.tau_sat_2 - self.strength_1) * self.sliprule.sum_slip(self.S, self.Q, self.H, self.L, self.T)

    hrates = np.array([hrate_1, hrate_2])

    hist = self.model.hist(self.S, self.Q, self.H, self.L, self.T, self.sliprule)

    self.assertTrue(np.allclose(hrates, hist))



