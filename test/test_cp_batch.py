#!/usr/bin/env python3

from neml import models, interpolate, elasticity, history
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, batch
from neml.math import rotations, tensors

import unittest
import numpy as np

class TestBatch(unittest.TestCase):
  def setUp(self):
    self.N = 10

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

    self.model = singlecrystal.SingleCrystalModel(self.kmodel, self.L, 
        miter = 120)

    self.orientations = rotations.random_orientations(self.N)

    self.D = np.array([0.01,-0.002,-0.003,0.012,-0.04,0.01])
    self.W = np.array([0.02,-0.02,0.03])

    self.T = 300.0
    self.dt = 2.0

  def test_batch_init(self):
    H1 = batch.init_history_batch(self.model, self.N)
    H2 = np.array([self.model.init_store() for i in range(self.N)])

    self.assertTrue(np.allclose(H1, H2))

  def test_batch_set_get(self):
    H = batch.init_history_batch(self.model, self.N)

    batch.set_orientation_passive_batch(self.model, H, self.orientations)

    nq = batch.get_orientation_passive_batch(self.model, H)

    for q1,q2 in zip(self.orientations,nq):
      self.assertTrue(np.allclose(q1.quat,q2.quat))

  def test_batch_no_threads(self):
    self.batch_run(1)

  def test_batch_threads(self):
    self.batch_run(2)

  def batch_run(self, nthreads):
    h_n = batch.init_history_batch(self.model, self.N)
    batch.set_orientation_passive_batch(self.model, h_n, self.orientations)

    d_np1 = np.array([self.D for i in range(self.N)])
    w_np1 = np.array([self.W for i in range(self.N)])

    T_np1 = np.array([self.T for i in range(self.N)])
    T_n = np.array([self.T for i in range(self.N)])

    d_n = np.zeros((self.N,6))
    w_n = np.zeros((self.N,3))

    s_n = np.zeros((self.N,6))
    u_n = np.zeros((self.N,))
    p_n = np.zeros((self.N,))

    s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = batch.evaluate_crystal_batch(
        self.model, d_np1, d_n, w_np1, w_n, T_np1, T_n, self.dt, 0.0,
        s_n, h_n, u_n, p_n, nthreads = nthreads)

    other_s_np1 = np.array([self.model.update_ld_inc(
      d_np1_i, d_n_i, w_np1_i, w_n_i, T_np1_i, T_n_i, self.dt, 0.0,
      s_n_i, h_n_i, u_n_i, p_n_i)[0] for d_np1_i, d_n_i, w_np1_i,
      w_n_i, T_np1_i, T_n_i, s_n_i, h_n_i, u_n_i, p_n_i in zip(
        d_np1, d_n, w_np1, w_n, T_np1, T_n, s_n, h_n, u_n, p_n)])

    self.assertTrue(np.allclose(s_np1, other_s_np1))

