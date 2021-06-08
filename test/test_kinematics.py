#!/usr/bin/env python3

from neml import history, interpolate, elasticity
from neml.math import tensors, rotations
from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, crystaldamage

import common
from common import differentiate
from nicediff import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonKinematics(object):
  def test_d_stress_rate_d_stress(self):
    nd = diff_symmetric_symmetric(lambda s: self.model.stress_rate(s, 
      self.d, self.w, self.Q, self.H, self.L, self.T, self.fixed), self.S)

    d = self.model.d_stress_rate_d_stress(self.S, self.d, self.w, self.Q, 
        self.H, self.L, self.T, self.fixed)

    print(d)

    self.assertTrue(d, nd)

  def test_d_stress_rate_d_d(self):
    def dfn(d):
      return self.model.stress_rate(self.S, d, self.w, self.Q, self.H, self.L, self.T, 
          self.fixed)

    nd = diff_symmetric_symmetric(dfn, self.d)
    d = self.model.d_stress_rate_d_d(self.S, self.d, self.w, self.Q, self.H, self.L,
        self.T, self.fixed)

    self.assertEqual(nd, d)

  def test_d_stress_rate_d_w(self):
    def dfn(w):
      return self.model.stress_rate(self.S, self.d, w, self.Q, self.H, self.L, self.T, self.fixed)

    nd = diff_symmetric_skew(dfn, self.w)
    d = self.model.d_stress_rate_d_w(self.S, self.d, self.w, self.Q, self.H,
        self.L, self.T, self.fixed)
   
    self.assertTrue(np.allclose(nd.data, d.data))

  def test_d_stress_rate_d_d_decouple(self):
    def dfn(d):
      fixed = self.model.decouple(self.S, d, self.w, self.Q, self.H, self.L, self.T,
          history.History())
      return self.model.stress_rate(self.S, d, self.w, self.Q, self.H, self.L, self.T, fixed)

    nd = diff_symmetric_symmetric(dfn, self.d) - self.model.d_stress_rate_d_d(self.S, self.d, 
        self.w, self.Q, self.H, self.L, self.T, self.fixed)
    d = self.model.d_stress_rate_d_d_decouple(self.S, self.d, self.w, self.Q, self.H, self.L,
        self.T, self.fixed)

    self.assertTrue(np.allclose(nd.data, d.data, atol = 1.0e-3))

  def test_d_stress_rate_d_w_decouple(self):
    def dfn(w):
      fixed = self.model.decouple(self.S, self.d, w, self.Q, self.H, self.L, self.T,
          history.History())
      return self.model.stress_rate(self.S, self.d, w, self.Q, self.H, self.L, self.T, fixed)

    nd = diff_symmetric_skew(dfn, self.w) - self.model.d_stress_rate_d_w(self.S, self.d,
        self.w, self.Q, self.H, self.L, self.T, self.fixed)
    d = self.model.d_stress_rate_d_w_decouple(self.S, self.d, self.w, self.Q, self.H,
        self.L, self.T, self.fixed)

    self.assertTrue(np.allclose(nd.data, d.data, rtol = 1.0e-3))

  def test_d_stress_rate_d_history(self):
    d = self.model.d_stress_rate_d_history(self.S, self.d, self.w, self.Q, self.H,
        self.L, self.T, self.fixed)
    nd = diff_symmetric_history(lambda h: self.model.stress_rate(self.S, self.d,
      self.w, self.Q, h, self.L, self.T, self.fixed), self.H)

    T = np.array(d).reshape(nd.shape, order = 'F')

    self.assertTrue(np.allclose(T, nd, rtol = 1.0e-3))

  def test_d_history_rate_d_stress(self):
    nd = diff_history_symmetric(lambda s: self.model.history_rate(s, self.d,
      self.w, self.Q, self.H, self.L, self.T, self.fixed), self.S)
    d = np.array(self.model.d_history_rate_d_stress(self.S, self.d, self.w,
        self.Q, self.H, self.L, self.T, self.fixed))

    self.assertTrue(np.allclose(nd, d.reshape(nd.shape)))

  def test_d_history_rate_d_d(self):
    nd = diff_history_symmetric(lambda d: self.model.history_rate(self.S, d,
      self.w, self.Q, self.H, self.L, self.T, self.fixed), self.d)
    d = np.array(self.model.d_history_rate_d_d(self.S, self.d, self.w, self.Q,
        self.H, self.L, self.T,self.fixed))

    self.assertTrue(np.allclose(nd, d.reshape(nd.shape)))

  def test_d_history_rate_d_w(self):
    nd = diff_history_skew(lambda w: self.model.history_rate(self.S, self.d,
      w, self.Q, self.H, self.L, self.T, self.fixed), self.w)
    d = np.array(self.model.d_history_rate_d_w(self.S, self.d, self.w, self.Q,
      self.H, self.L, self.T, self.fixed))

    self.assertTrue(np.allclose(nd, d.reshape(nd.shape)))

  def test_d_history_rate_d_d_decouple(self):
    def dfn(d):
      fixed = self.model.decouple(self.S, d, self.w, self.Q, self.H, self.L, self.T,
          history.History())
      return self.model.history_rate(self.S, d, self.w, self.Q, self.H, self.L, self.T, fixed)
    
    nd1 = diff_history_symmetric(dfn, self.d)

    nd = nd1 - np.array(self.model.d_history_rate_d_d(self.S, self.d, self.w,
        self.Q, self.H, self.L, self.T, self.fixed)).reshape(nd1.shape, order = 'F')
    d = np.array(self.model.d_history_rate_d_d_decouple(self.S, self.d, self.w, self.Q,
        self.H, self.L, self.T, self.fixed))

    self.assertTrue(np.allclose(nd, d.reshape(nd.shape, order = 'F')))
  
  def test_d_history_rate_d_w_decouple(self):
    def dfn(w):
      fixed = self.model.decouple(self.S, self.d, w, self.Q, self.H, self.L, self.T,
          history.History())
      return self.model.history_rate(self.S, self.d, w, self.Q, self.H, self.L, self.T, fixed)
  
    nd1 = diff_history_skew(dfn, self.w)

    nd = nd1 - np.array(self.model.d_history_rate_d_w(self.S, self.d, self.w,
        self.Q, self.H, self.L, self.T, self.fixed)).reshape(nd1.shape, order = 'F')
    d = np.array(self.model.d_history_rate_d_w_decouple(self.S, self.d, self.w, self.Q,
        self.H, self.L, self.T, self.fixed))

    self.assertTrue(np.allclose(nd, d.reshape(nd.shape, order = 'F')))

  def test_d_history_rate_d_history(self):
    nd = diff_history_history(lambda h: self.model.history_rate(self.S, self.d,
      self.w, self.Q, h, self.L, self.T, self.fixed), self.H)


    d = np.array(self.model.d_history_rate_d_history(self.S, self.d, self.w,
      self.Q, self.H, self.L, self.T, self.fixed))

    T = d.reshape(nd.shape)
    
    self.assertTrue(np.allclose(nd, T))

class TestStandardKinematics(unittest.TestCase, CommonKinematics):
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

    self.imodel = inelasticity.AsaroInelasticity(self.slipmodel)

    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))

    self.T = 300.0

    self.mu = 29000.0
    self.E = 120000.0
    self.nu = 0.3

    self.emodel = elasticity.CubicLinearElasticModel(self.E, 
        self.nu, self.mu, "moduli")

    self.dn = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.dn = 0.5*(self.dn + self.dn.T)
    self.d = tensors.Symmetric(self.dn)

    self.wn = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.wn = 0.5 * (self.wn - self.wn.T)
    self.w = tensors.Skew(self.wn)

    self.model = kinematics.StandardKinematicModel(self.emodel, self.imodel)
    
    self.fspin = self.model.spin(self.S, self.d, self.w, self.Q, self.H,
        self.L, self.T, history.History())

    self.fixed = self.model.decouple(self.S, self.d, self.w, self.Q, self.H, self.L, self.T,
        history.History())

  def test_setup_history(self):
    H1 = history.History()
    self.model.populate_history(H1)
    self.model.init_history(H1)

    H2 = history.History()
    self.imodel.populate_history(H2)
    self.imodel.init_history(H2)

    self.assertTrue(np.allclose(np.array(H1), np.array(H2)))

  def test_stress_rate(self):
    Cfull = ms2ts(self.emodel.C_tensor(self.T, self.Q).data.reshape((6,6)))
    Sfull = ms2ts(self.emodel.S_tensor(self.T, self.Q).data.reshape((6,6)))
    d = usym(self.d.data)
    w = uskew(self.w.data)
    Ofull = uskew(self.fspin.data)
    
    dp = usym(self.imodel.d_p(self.S, self.Q, self.H, self.L, self.T, self.fixed).data)
    wp = uskew(self.imodel.w_p(self.S, self.Q, self.H, self.L, self.T, self.fixed).data)

    O = wp + Ofull

    stress = usym(self.S.data)
    
    e = np.einsum('ijkl,kl', Sfull, stress)

    sdot1 = tensors.Symmetric(
        np.einsum('ijkl,kl', Cfull, d - dp - np.dot(e, O) + np.dot(O, e)))
    
    sdot2 = self.model.stress_rate(self.S, self.d, self.w, self.Q, self.H,
        self.L, self.T, self.fixed)

    self.assertEqual(sdot1, sdot2)

  def test_hist_rate(self):
    H1 = self.model.history_rate(self.S, self.d, self.w,
        self.Q, self.H, self.L, self.T, self.fixed)
    H2 = self.imodel.history_rate(self.S, self.Q, self.H, self.L, self.T, self.fixed)

    self.assertTrue(np.allclose(np.array(H1), np.array(H2)))

  def test_spin_rate(self):
    Cfull = ms2ts(self.emodel.C_tensor(self.T, self.Q).data.reshape((6,6)))
    Sfull = ms2ts(self.emodel.S_tensor(self.T, self.Q).data.reshape((6,6)))
    d = usym(self.d.data)
    w = uskew(self.w.data)
    Ofull = uskew(self.fspin.data)
    
    dp = usym(self.imodel.d_p(self.S, self.Q, self.H, self.L, self.T, self.fixed).data)
    wp = uskew(self.imodel.w_p(self.S, self.Q, self.H, self.L, self.T, self.fixed).data)

    O = wp + Ofull

    stress = usym(self.S.data)
    
    e = np.einsum('ijkl,kl', Sfull, stress)

    spin1 = tensors.Skew(
        w - wp - np.dot(e, dp) + np.dot(dp, e))
    spin2 = self.model.spin(self.S, self.d, self.w, self.Q, self.H, 
        self.L, self.T, self.fixed)

    self.assertTrue(spin1, spin2)

class TestStandardKinematicsComplicated(unittest.TestCase, CommonKinematics):
  def setUp(self):
    self.strength_0 = 35.0
    self.H = history.History()
    self.H.add_scalar("strength0")
    self.H.set_scalar("strength0", self.strength_0)

    self.strength_1 = 25.0
    self.H.add_scalar("strength1")
    self.H.set_scalar("strength1", self.strength_1)

    self.tau0_0 = 10.0
    self.tau_sat_0 = 50.0
    self.b_0 = 2.5

    self.tau0_1 = 5.0
    self.tau_sat_1 = 25.0
    self.b_1 = 1.0

    self.strengthmodel = slipharden.SumSlipSingleStrengthHardening(
        [
          slipharden.VoceSlipHardening(self.tau_sat_0, self.b_0, self.tau0_0),
          slipharden.VoceSlipHardening(self.tau_sat_1, self.b_1, self.tau0_1)
          ])
    
    self.g0 = 1.0
    self.n = 3.0
    self.slipmodel = sliprules.PowerLawSlipRule(self.strengthmodel, self.g0, self.n)

    self.imodel = inelasticity.AsaroInelasticity(self.slipmodel)

    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))

    self.T = 300.0

    self.mu = 29000.0
    self.E = 120000.0
    self.nu = 0.3

    self.emodel = elasticity.CubicLinearElasticModel(self.E, 
        self.nu, self.mu, "moduli")

    self.dn = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.dn = 0.5*(self.dn + self.dn.T)
    self.d = tensors.Symmetric(self.dn)

    self.wn = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.wn = 0.5 * (self.wn - self.wn.T)
    self.w = tensors.Skew(self.wn)

    self.model = kinematics.StandardKinematicModel(self.emodel, self.imodel)
    
    self.fspin = self.model.spin(self.S, self.d, self.w, self.Q, self.H,
        self.L, self.T, history.History())

    self.fixed = self.model.decouple(self.S, self.d, self.w, self.Q, self.H, self.L, self.T,
        history.History())

class TestDamagedStandardKinematics(unittest.TestCase, CommonKinematics):
  def setUp(self):
    self.strength = 35.0
    self.H = history.History()
    self.H.add_scalar("strength")
    self.H.set_scalar("strength", self.strength)
    self.H.add_scalar("whatever")
    self.H.set_scalar("whatever", 0.5)

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
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))

    self.T = 300.0

    self.mu = 29000.0
    self.E = 120000.0
    self.nu = 0.3

    self.emodel = elasticity.CubicLinearElasticModel(self.E, 
        self.nu, self.mu, "moduli")

    self.dn = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.dn = 0.5*(self.dn + self.dn.T)
    self.d = tensors.Symmetric(self.dn)

    self.wn = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.wn = 0.5 * (self.wn - self.wn.T)
    self.w = tensors.Skew(self.wn)

    self.dmodel = crystaldamage.NilDamageModel()

    self.model = kinematics.DamagedStandardKinematicModel(self.emodel, self.imodel,
        self.dmodel)
    
    self.fspin = self.model.spin(self.S, self.d, self.w, self.Q, self.H,
        self.L, self.T, history.History())

    self.fixed = self.model.decouple(self.S, self.d, self.w, self.Q, self.H, self.L, self.T,
        history.History())

class TestDamagedStandardKinematicsComplex(unittest.TestCase, CommonKinematics):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])

    self.nslip = self.L.ntotal

    self.strength = 35.0
    self.c = 10.0
    self.beta  = 2.0

    self.H = history.History()

    for i in range(12):
      self.H.add_scalar("strength"+str(i))
      self.H.set_scalar("strength"+str(i), self.strength)
    
    for j in range(4):
      self.H.add_scalar("slip_damage_"+str(j))
      self.H.set_scalar("slip_damage_"+str(j), self.c*0.4) 

    self.static = 20.0
    self.s0 = [self.static]*self.nslip

    self.k = 1000.0
    self.sat = 40.0
    self.m = 1.5

    self.strengthmodel = slipharden.VocePerSystemHardening(
        self.s0, [self.k]*self.nslip, [self.sat]*self.nslip,
        [self.m]*self.nslip)

    self.g0 = 1.0
    self.n = 3.0
    self.slipmodel = sliprules.PowerLawSlipRule(self.strengthmodel, self.g0, self.n)

    self.imodel = inelasticity.AsaroInelasticity(self.slipmodel)


    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))

    self.T = 300.0

    self.mu = 29000.0
    self.E = 120000.0
    self.nu = 0.3

    self.emodel = elasticity.CubicLinearElasticModel(self.E, 
        self.nu, self.mu, "moduli")

    self.dn = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.dn = 0.5*(self.dn + self.dn.T)
    self.d = tensors.Symmetric(self.dn)

    self.wn = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.wn = 0.5 * (self.wn - self.wn.T)
    self.w = tensors.Skew(self.wn)

    self.dmodel = crystaldamage.WorkPlaneDamage()
    self.nfunc  = crystaldamage.SigmoidTransformation(self.c, self.beta)
    self.sfunc  = crystaldamage.SigmoidTransformation(self.c, self.beta)

    self.dmodel = crystaldamage.PlanarDamageModel(self.dmodel, self.nfunc, self.sfunc,
        self.L)

    self.model = kinematics.DamagedStandardKinematicModel(self.emodel, self.imodel,
        self.dmodel)
    
    self.fspin = self.model.spin(self.S, self.d, self.w, self.Q, self.H,
        self.L, self.T, history.History())

    self.fixed = self.model.decouple(self.S, self.d, self.w, self.Q, self.H, self.L, self.T,
        history.History())

     
