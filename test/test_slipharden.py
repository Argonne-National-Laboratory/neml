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
    d = np.array(self.model.d_hist_d_s(self.S, self.Q, self.H, self.L, self.T, self.sliprule, self.fixed))
    nd = diff_history_symmetric(lambda s: self.model.hist(s, self.Q, self.H, self.L, self.T,
      self.sliprule, self.fixed), self.S)
    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

  def test_d_hist_d_hist(self):
    d = np.array(self.model.d_hist_d_h(self.S, self.Q, self.H, self.L, self.T, self.sliprule, 
      self.fixed))
    nd = diff_history_history(lambda h: self.model.hist(self.S, self.Q, h, self.L, self.T,
      self.sliprule, self.fixed), self.H)
    
    print(d)
    print(nd)
    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

  def test_d_hist_to_tau_d_hist(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        nd = diff_history_scalar(lambda h: self.model.hist_to_tau(g, i, h, self.T, self.fixed), self.H)
        d = self.model.d_hist_to_tau(g, i, self.H, self.T, self.fixed)
        self.assertTrue(np.allclose(np.array(nd), np.array(d)))

class CommonSlipSingleHardening():
  def test_d_hist_map(self):
    nd = diff_history_scalar(lambda h: self.model.hist_map(h, self.T, 
      self.fixed), self.H)
    H = self.model.d_hist_map(self.H, self.T, self.fixed)
    self.assertTrue(np.allclose(np.array(H), np.array(nd)))

  def test_hist_to_tau(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        self.assertTrue(np.isclose(self.strength + self.static + self.nye_part, 
          self.model.hist_to_tau(g, i, self.H, self.T, self.fixed)))

class CommonSlipSingleStrengthHardening():
  def test_hist_map(self):
    self.assertTrue(np.isclose(self.model.hist_map(self.H, self.T,
      self.fixed), self.strength + self.static + self.nye_part))

  def test_hist(self):
    h = self.model.hist(self.S, self.Q, self.H, self.L, self.T, self.sliprule, self.fixed)
    self.assertTrue(np.isclose(h.get_scalar("strength"), self.model.hist_rate(
      self.S, self.Q, self.H, self.L, self.T, self.sliprule, self.fixed)))

  def test_d_hist_rate_d_stress(self):
    dfn = lambda s: self.model.hist_rate(s, self.Q, self.H,
        self.L, self.T, self.sliprule, self.fixed)
    nd = diff_scalar_symmetric(dfn, self.S)
    d = self.model.d_hist_rate_d_stress(self.S, self.Q, self.H, self.L,
        self.T, self.sliprule, self.fixed)
    self.assertEqual(d,nd)

  def test_d_hist_rate_d_hist(self):
    dfn = lambda h: self.model.hist_rate(self.S, self.Q, h, 
        self.L, self.T, self.sliprule, self.fixed)
    nd = diff_history_scalar(dfn, self.H)
    
    d = self.model.d_hist_rate_d_hist(self.S, self.Q, self.H, self.L,
        self.T, self.sliprule, self.fixed)
    
    self.assertTrue(np.isclose(d.get_scalar("strength"), nd.get_scalar("strength")))

  def test_nye_parts_match(self):
    self.assertTrue(np.isclose(self.model.nye_part(self.nye, self.T), 
      self.model.nye_contribution(self.fixed, self.T)))

  def test_no_nye_works(self):
    self.assertTrue(np.isclose(self.model.nye_contribution(history.History(), self.T),
      0.0))

class CommonPlasticSlipHardening():
  def sum_slip(self):
    return sum(np.abs(self.sliprule.slip(g, i, self.S, self.Q, self.H, self.L, 
      self.T, self.fixed)) for g in range(self.L.ngroup) for i in range(self.L.nslip(g)))

  def test_hist_rate(self):
    ss = self.sum_slip()

    self.assertTrue(np.isclose(self.model.hist_rate(self.S, self.Q, self.H,
      self.L, self.T, self.sliprule, self.fixed), self.model.hist_factor(self.H.get_scalar("strength"), 
        self.L, self.T, self.fixed) * ss))

  def test_d_factor(self):
    nd = differentiate(lambda s: self.model.hist_factor(s, self.L, self.T, self.fixed), 
        self.strength)
    d = self.model.d_hist_factor(self.strength, self.L, self.T, self.fixed)

    self.assertTrue(np.isclose(nd, d))
    
class TestVoceHardening(unittest.TestCase, CommonPlasticSlipHardening, 
    CommonSlipSingleStrengthHardening, CommonSlipSingleHardening, CommonSlipHardening):
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

    self.static = self.tau0

    self.model = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0)
    
    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)

    self.fixed = history.History()

    self.nye = tensors.RankTwo([[1.1,1.2,1.3],[2.1,2.2,2.3],[3.1,3.2,3.3]])

    self.nye_part = 0.0

  def test_initialize_hist(self):
    H = history.History()
    self.model.populate_history(H)
    self.model.init_history(H)
    self.assertTrue(np.isclose(H.get_scalar('strength'), 0.0))

  def test_static_strength(self):
    self.assertTrue(np.isclose(self.model.static_strength(self.T), self.tau0))
  
  def test_factor(self):
    self.assertTrue(np.isclose(
      self.model.hist_factor(self.strength, self.L, self.T, self.fixed),
      self.b * (self.tau_sat - self.strength)))

  def test_nye_factor(self):
    self.assertTrue(np.isclose(
      self.model.nye_part(self.nye, self.T), 0.0))

  def test_use_nye(self):
    self.assertFalse(self.model.use_nye)

class TestNyeVoceHardening(unittest.TestCase, CommonPlasticSlipHardening, 
    CommonSlipSingleStrengthHardening, CommonSlipSingleHardening, CommonSlipHardening):
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

    self.static = self.tau0

    self.k = 1.2

    self.model = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0, 
        k = self.k)
    
    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)

    self.nye = tensors.RankTwo([[1.1,1.2,1.3],[2.1,2.2,2.3],[3.1,3.2,3.3]])
    self.nye_part = self.k * np.sqrt(la.norm(self.nye.data.reshape((3,3)), ord = 'fro'))

    self.fixed = history.History()
    self.fixed.add_ranktwo("nye")
    self.fixed.set_ranktwo("nye", self.nye)
    
  def test_initialize_hist(self):
    H = history.History()
    self.model.populate_history(H)
    self.model.init_history(H)
    self.assertTrue(np.isclose(H.get_scalar('strength'), 0.0))

  def test_static_strength(self):
    self.assertTrue(np.isclose(self.model.static_strength(self.T), self.tau0))
  
  def test_factor(self):
    self.assertTrue(np.isclose(
      self.model.hist_factor(self.strength, self.L, self.T, self.fixed),
      self.b * (self.tau_sat - self.strength)))

  def test_nye_factor(self):
    self.assertTrue(np.isclose(
      self.model.nye_part(self.nye, self.T), self.nye_part))

  def test_use_nye(self):
    self.assertTrue(self.model.use_nye)

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

    self.fixed = history.History()

    self.nye = tensors.RankTwo([[1.1,1.2,1.3],[2.1,2.2,2.3],[3.1,3.2,3.3]])

    self.nye_part = 0.0

  def test_initialize_hist(self):
    H = history.History()
    self.model.populate_history(H)
    self.model.init_history(H)
    self.assertTrue(np.isclose(H.get_scalar('strength0'), 0.0))
    self.assertTrue(np.isclose(H.get_scalar('strength1'), 0.0))

  def test_hist_def(self):
    hrate_1 = self.b_1 * (self.tau_sat_1 - self.strength_0) * self.sliprule.sum_slip(self.S, self.Q, 
        self.H, self.L, self.T, self.fixed)
    hrate_2 = self.b_2 * (self.tau_sat_2 - self.strength_1) * self.sliprule.sum_slip(self.S, self.Q,
        self.H, self.L, self.T, self.fixed)

    hrates = np.array([hrate_1, hrate_2])

    hist = self.model.hist(self.S, self.Q, self.H, self.L, self.T, self.sliprule, self.fixed)

    self.assertTrue(np.allclose(hrates, hist))

  def test_use_nye(self):
    self.assertFalse(self.model.use_nye)

class TestLinearSlipHardening(unittest.TestCase, CommonPlasticSlipHardening, 
    CommonSlipSingleStrengthHardening, CommonSlipSingleHardening, CommonSlipHardening):
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
    self.k1 = 2.0
    self.k2 = 3.0

    self.static = self.tau0

    self.model = slipharden.LinearSlipHardening(self.tau0, self.k1, self.k2)
    
    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)

    self.nye = tensors.RankTwo([[1.1,1.2,1.3],[2.1,2.2,2.3],[3.1,3.2,3.3]])
    self.nye_part = self.k2 * la.norm(self.nye.data.reshape((3,3)), ord = 'fro')

    self.fixed = history.History()
    self.fixed.add_ranktwo("nye")
    self.fixed.set_ranktwo("nye", self.nye)
    
  def test_initialize_hist(self):
    H = history.History()
    self.model.populate_history(H)
    self.model.init_history(H)
    self.assertTrue(np.isclose(H.get_scalar('strength'), 0.0))

  def test_static_strength(self):
    self.assertTrue(np.isclose(self.model.static_strength(self.T), self.tau0))
  
  def test_factor(self):
    self.assertTrue(np.isclose(
      self.model.hist_factor(self.strength, self.L, self.T, self.fixed),
      self.k1))

  def test_nye_factor(self):
    self.assertTrue(np.isclose(
      self.model.nye_part(self.nye, self.T), self.nye_part))

  def test_use_nye(self):
    self.assertTrue(self.model.use_nye)

