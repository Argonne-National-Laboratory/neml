#!/usr/bin/env python3

from neml import history, interpolate
from neml.math import tensors, rotations
from neml.cp import crystallography, slipharden, sliprules

from common import differentiate
from nicediff import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonSlipRule(object):
  def test_d_slip_d_stress(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        d = self.model.d_slip_d_s(g, i, self.S, self.Q, self.H, self.L, self.T)
        nd = diff_scalar_symmetric(lambda s: self.model.slip(g, i, s, self.Q, self.H, 
          self.L, self.T), self.S)
        self.assertEqual(d, nd)

  def test_d_slip_d_hist(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        d = np.array(self.model.d_slip_d_h(g, i, self.S, self.Q, self.H, self.L, self.T))
        nd = np.array(diff_history_scalar(lambda h: self.model.slip(g, i, self.S, self.Q, h,
          self.L, self.T), self.H))
        self.assertTrue(np.allclose(nd.reshape(d.shape), d))

  def test_d_hist_rate_d_stress(self):
    d = np.array(self.model.d_hist_rate_d_stress(self.S, self.Q, self.H, self.L, self.T))
    nd = diff_history_symmetric(lambda s: self.model.hist_rate(s, self.Q, self.H, self.L,
      self.T), self.S)
    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

  def test_d_hist_rate_d_hist(self):
    d = np.array(self.model.d_hist_rate_d_hist(self.S, self.Q, self.H, self.L, self.T))
    nd = diff_history_history(lambda h: self.model.hist_rate(self.S, self.Q, h, self.L,
      self.T), self.H)
    self.assertTrue(np.allclose(nd.reshape(d.shape), d))

class CommonSlipMultiStrengthSlipRule(object):
  def test_setup_history(self):
    model_hist = history.History()
    self.model.populate_history(model_hist)
    self.model.init_history(model_hist)
    
    manual_hist = history.History()
    for strength in self.strengths:
      Hs = history.History()
      strength.populate_history(Hs)
      strength.init_history(Hs)
      manual_hist.add_union(Hs)

    self.assertTrue(np.allclose(np.array(model_hist),
      np.array(manual_hist)))

  def test_history_rate(self):
    reference = self.model.hist_rate(self.S, self.Q, self.H, self.L, self.T)
    combined = history.History()
    for strength in self.strengths:
      combined.add_union(strength.hist(self.S, self.Q, self.H, self.L, self.T, self.model))

    self.assertTrue(np.allclose(reference, combined))

  def test_slip_rate(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        rs = self.L.shear(g, i, self.Q, self.S)
        print(self.model.slip(g,i,self.S, self.Q, self.H, self.L, self.T))
        print(self.model.sslip(g, i, rs, self.strength_values, self.T))
        print(self.strength_values)
        self.assertTrue(np.isclose(self.model.slip(g, i, self.S, self.Q, self.H, self.L, self.T),
          self.model.sslip(g, i, rs, self.strength_values, self.T)))

  def test_slip_tau(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        rs = self.L.shear(g, i, self.Q, self.S)
        nd = differentiate(lambda t: self.model.sslip(g, i, t, self.strength_values, self.T),
            rs)
        d = self.model.d_sslip_dtau(g, i, rs, self.strength_values, self.T)

        self.assertTrue(np.isclose(d, nd))

  def test_slip_strength(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        rs = self.L.shear(g, i, self.Q, self.S)

        nd = differentiate(lambda s: self.model.sslip(g, i, rs, s, self.T), np.array(self.strength_values))
        d = self.model.d_sslip_dstrength(g, i, rs, self.strength_values, self.T)

        self.assertTrue(np.allclose(nd,d))

class TestKinematicPowerLawSlip(unittest.TestCase, CommonSlipMultiStrengthSlipRule):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))
    self.strength0 = 15.0
    self.strength1 = 20.0
    self.H = history.History()
    self.H.add_scalar("strength0")
    self.H.add_scalar("strength1")
    self.H.set_scalar("strength0", self.strength0)
    self.H.set_scalar("strength1", self.strength1)

    self.T = 300.0

    self.tau0_1 = 10.0
    self.tau_sat_1 = 50.0
    self.b_1 = 2.5

    self.tau0_2 = 15.0
    self.tau_sat_2 = 55.0
    self.b_2 = 3.0

    self.backstrength = slipharden.VoceSlipHardening(self.tau_sat_1, self.b_1, self.tau0_1)
    self.understrength = slipharden.VoceSlipHardening(self.tau_sat_2, self.b_2, self.tau0_2)
    self.strengths = [self.backstrength, self.understrength]

    self.strength_values = [self.strength0 + self.tau0_1, self.strength1 + self.tau0_2]
    
    self.g0 = 1.0
    self.n = 3.0
    self.model = sliprules.KinematicPowerLawSlipRule(self.backstrength, self.understrength, self.g0, self.n)

    self.tau = 45.0

  def test_sslip(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        inv = (self.tau - self.strength0 - self.tau0_1) / (self.strength1 + self.tau0_2) 

        expected = self.g0 * np.abs(inv) ** (self.n-1.0) * inv
        
        actual = self.model.sslip(g, i, self.tau, self.strength_values, self.T)

        self.assertTrue(np.isclose(expected, actual))

class CommonSlipStrengthSlipRule(CommonSlipMultiStrengthSlipRule):
  def test_init_hist(self):
    H1 = history.History()
    self.model.populate_history(H1)
    self.model.init_history(H1)

    H2 = history.History()
    self.strengthmodel.populate_history(H2)
    self.strengthmodel.init_history(H2)

    self.assertTrue(np.allclose(np.array(H1),
      np.array(H2)))

  def test_slip(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        rs = self.L.shear(g, i, self.Q, self.S)
        strength = self.strength + self.static
        self.assertTrue(np.isclose(self.model.slip(g, i, self.S, self.Q, self.H, self.L, self.T),
          self.model.scalar_sslip(g, i, rs, strength, self.T)))

  def test_d_hist_rate(self):
    self.assertTrue(np.allclose(
      np.array(self.model.hist_rate(self.S, self.Q, self.H, self.L, self.T)),
      np.array(self.strengthmodel.hist(self.S, self.Q, self.H, self.L, self.T, self.model))))

  def test_d_sslip_d_tau(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        nd = differentiate(lambda t: self.model.scalar_sslip(g, i, t, self.strength, self.T),
            self.tau)
        d = self.model.scalar_d_sslip_dtau(g, i, self.tau, self.strength, self.T)
        self.assertTrue(np.isclose(nd,d))

  def test_d_sslip_d_strength(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        nd = differentiate(lambda s: self.model.scalar_sslip(g, i, self.tau, s, self.T), self.strength)
        d = self.model.scalar_d_sslip_dstrength(g, i, self.tau, self.strength, self.T)
        self.assertTrue(np.isclose(nd, d))

class TestPowerLawSlip(unittest.TestCase, CommonSlipStrengthSlipRule, CommonSlipRule):
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

    self.strengthmodel = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0)
    self.strengths = [self.strengthmodel]

    self.static = self.tau0

    self.strength_values = [self.strength + self.static]
    
    self.g0 = 1.0
    self.n = 3.0
    self.model = sliprules.PowerLawSlipRule(self.strengthmodel, self.g0, self.n)

    self.tau = 33.0

  def test_scalar_rate(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        self.assertTrue(np.isclose(self.model.scalar_sslip(g, i, self.tau, self.strength, self.T),
          self.g0 * np.abs(self.tau/self.strength)**(self.n-1.0) * self.tau/self.strength))


class TestBiVoceSlip(unittest.TestCase, CommonSlipStrengthSlipRule, CommonSlipRule):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))
    self.strength_1 = 35.0
    self.strength_2 = 25.0
    self.strength = self.strength_1 + self.strength_2
    self.H = history.History()
    self.H.add_scalar("strength0")
    self.H.set_scalar("strength0", self.strength_1)
    self.H.add_scalar("strength1")
    self.H.set_scalar("strength1", self.strength_2)

    self.T = 300.0

    self.tau0 = 10.0
    self.tau_sat = 50.0
    self.b = 2.5

    self.strengthmodel = slipharden.SumSlipSingleStrengthHardening(
        [slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0),
          slipharden.VoceSlipHardening(self.tau_sat/2, self.b/2, self.tau0/2)])
    self.static = self.tau0 + self.tau0 / 2

    self.strengths = [self.strengthmodel]
    self.strength_values = [self.strength + self.static]
    
    self.g0 = 1.0
    self.n = 3.0
    self.model = sliprules.PowerLawSlipRule(self.strengthmodel, self.g0, self.n)

    self.tau = 33.0
