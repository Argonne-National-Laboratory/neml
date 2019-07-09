#!/usr/bin/env python3

from neml import models, interpolate, elasticity, history
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal
from neml.math import rotations, tensors

import unittest
import numpy as np

import common

class CommonSolver(object):
  pass

class TestSingleCrystal(unittest.TestCase, CommonSolver):
  def setUp(self):
    self.tau0 = 10.0
    self.tau_sat = 50.0
    self.b = 2.5

    self.strengthmodel = slipharden.VoceSlipHardening(self.tau_sat, self.b, self.tau0)
    
    self.g0 = 1.0
    self.n = 3.0
    self.slipmodel = sliprules.PowerLawSlipRule(self.strengthmodel, self.g0, self.n)

    self.imodel = inelasticity.AsaroInelasticity(self.slipmodel)

    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")

    self.mu = 29000.0
    self.E = 120000.0
    self.nu = 0.3

    self.emodel = elasticity.CubicLinearElasticModel(self.E, 
        self.nu, self.mu, "moduli")

    self.kmodel = kinematics.StandardKinematicModel(self.emodel, self.imodel)

    self.model = singlecrystal.SingleCrystalModel(self.kmodel, self.L, initial_rotation = self.Q)

    self.T = 300.0
    self.stress_n = np.array([120.0,-60.0,170.0,35.0,80.0,-90.0])
    self.stress_np1 = np.array([15.0,-40.0,120.0,70.0,-10.0,-50.0])

    self.d = np.array([0.1,-0.2,0.25,0.11,-0.05,0.075])
    self.w = np.array([0.1,0.2,-0.2])

    self.strength_n = 25.0
    self.strength_np1 = 30.0

    self.S_np1 = tensors.Symmetric(common.usym(self.stress_np1))
    self.S_n = tensors.Symmetric(common.usym(self.stress_n))

    self.D = tensors.Symmetric(common.usym(self.d))
    self.W = tensors.Skew(common.uskew(self.w))

    self.strength_n = 25.0
    self.strength_np1 = 30.0

    self.H_n = history.History()
    self.H_n.add_scalar("strength")
    self.H_n.set_scalar("strength", self.strength_n)

    self.H_np1 = history.History()
    self.H_np1.add_scalar("strength")
    self.H_np1.set_scalar("strength", self.strength_np1)

    self.dt = 2.0

    self.ts = singlecrystal.SCTrialState(self.D, self.W, self.S_n, self.H_n, self.Q, self.L, self.T, self.dt)

    self.x = np.zeros((self.model.nparams,))
    self.x[:6] = self.stress_np1
    self.x[6] = self.strength_np1

  def test_nhist(self):
    self.assertEqual(self.model.nstore, 5)

  def test_init_hist(self):
    h = self.model.init_store()
    
    self.assertTrue(np.allclose(h[:4], self.Q.quat))
    self.assertTrue(np.isclose(h[4],0.0))

  def test_residual(self):
    R, J = self.model.RJ(self.x, self.ts)

    Rtrue = np.zeros((self.model.nparams,))
    srate = self.kmodel.stress_rate(self.S_np1, self.D, self.W, self.Q, self.H_np1, self.L, self.T)

    hrate = np.array(self.kmodel.history_rate(self.S_np1, self.D, self.W, self.Q, self.H_np1, self.L, self.T))

    Rtrue[:6] = self.stress_np1 - self.stress_n - srate.data * self.dt
    Rtrue[6:] = np.array(self.H_np1) - np.array(self.H_n) - hrate * self.dt

    self.assertTrue(np.allclose(Rtrue, R))
