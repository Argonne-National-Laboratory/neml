#!/usr/bin/env python3

from neml import history, interpolate
from neml.math import tensors, rotations, matrix, projections
from neml.cp import crystallography, sliprules, crystaldamage, slipharden

from common import differentiate
from nicediff import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonCrystalDamageModel():
  def test_projection_stress(self):
    dfn = lambda s: self.model.projection(s, self.huse, self.Q,
        self.L, self.sliprule, self.T)
    R1 = diff_symsym_sym(dfn, self.S)
    R2 = self.model.d_projection_d_stress(self.S, self.huse, self.Q,
        self.L, self.sliprule, self.T)

    self.assertEqual(R1,R2)

  def test_projection_history(self):
    dfn = lambda h: self.model.projection(self.S, h, self.Q,
        self.L, self.sliprule, self.T)
    res = diff_symsym_history(dfn, self.hdmg)
    act = np.array(self.model.d_projection_d_history(self.S, self.hdmg, self.Q,
        self.L, self.sliprule, self.T))
    
    shape = res.shape[::-1]

    self.assertTrue(np.allclose(res, act.reshape(shape).T)) # Transpose?

  def test_history_stress(self):
    dd = diff_history_symmetric(
        lambda s: self.model.damage_rate(s, self.huse, self.Q,
        self.L, self.sliprule, self.T, self.fixed), self.S)
    d = np.array(self.model.d_damage_d_stress(self.S, self.huse, self.Q,
      self.L, self.sliprule, self.T, self.fixed))

    self.assertTrue(np.allclose(dd.flatten(), d))

  def test_history_history(self):
    nbase = len(np.array(self.hbase))
    ndmg = len(np.array(self.hdmg))
    ntotal = nbase + ndmg

    d = np.array(self.model.d_damage_d_history(self.S, self.huse, self.Q,
      self.L, self.sliprule, self.T, self.fixed)).reshape(ndmg,ntotal).flatten()
    nd = diff_history_history(lambda h: self.model.damage_rate(
      self.S, h, self.Q, self.L, self.sliprule, self.T, self.fixed), self.huse)

    self.assertTrue(np.allclose(d, nd.flatten()))

class TestPlanarDamageModel(unittest.TestCase, CommonCrystalDamageModel):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])

    self.nslip = self.L.ntotal
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))*2
    
    self.static = 20.0
    self.s0 = [self.static]*self.nslip

    self.k = 1000.0
    self.sat = 40.0
    self.m = 1.5

    self.hmodel = slipharden.VocePerSystemHardening(
        self.s0, [self.k]*self.nslip, [self.sat]*self.nslip,
        [self.m]*self.nslip)

    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.hmodel, self.g0, self.n)

    self.T = 300.0

    self.c = 10.0
    self.beta  = 2.0

    self.dmodel = crystaldamage.WorkPlaneDamage()
    self.nfunc  = crystaldamage.SigmoidTransformation(self.c, self.beta)
    self.sfunc  = crystaldamage.SigmoidTransformation(self.c, self.beta)

    self.model = crystaldamage.PlanarDamageModel(self.dmodel, self.sfunc,
        self.nfunc, self.L)

    self.huse = history.History()
    self.hmodel.populate_history(self.huse)
    self.model.populate_history(self.huse)

    for i in range(12):
      self.huse.set_scalar("strength"+str(i), 2.0)
    
    for j in range(4):
      self.huse.set_scalar("slip_damage_"+str(j), self.c*0.4) 

    self.hbase = self.huse.subset(["strength"+str(i) for i in range(12)])
    self.hdmg  = self.huse.subset(["slip_damage_"+str(i) for i in range(4)])

    self.fixed = history.History()

  def test_projection(self):
    should = tensors.SymSymR4.id()
    
    for i in range(self.L.nplanes):
      n = self.Q.apply(self.L.unique_planes[i])
      
      Pn = projections.normal_projection_ss(n)
      Ps = projections.shear_projection_ss(n)

      ns = self.S.dot(n).dot(n)
      d = self.huse.get_scalar("slip_damage_"+str(i))

      N = self.nfunc.map(d, ns)
      S = self.sfunc.map(d, ns)
      
      should = should.dot(tensors.SymSymR4.id() - N*Pn - S*Ps)

    actual = self.model.projection(self.S, self.huse, self.Q,
        self.L, self.sliprule, self.T)
    self.assertTrue(np.allclose(actual.data,should.data))

class TestNilDamageModel(unittest.TestCase, CommonCrystalDamageModel):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])

    self.nslip = self.L.ntotal
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))
    
    self.static = 20.0
    self.s0 = [self.static]*self.nslip

    self.k = 1000.0
    self.sat = 40.0
    self.m = 1.5

    self.hmodel = slipharden.VocePerSystemHardening(
        self.s0, [self.k]*self.nslip, [self.sat]*self.nslip,
        [self.m]*self.nslip)

    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.hmodel, self.g0, self.n)

    self.T = 300.0

    self.model = crystaldamage.NilDamageModel()

    self.huse = history.History()
    self.hmodel.populate_history(self.huse)
    self.model.populate_history(self.huse)

    for i in range(12):
      self.huse.set_scalar("strength"+str(i), 25.0)

    self.huse.set_scalar("whatever", 0.5)

    self.hbase = self.huse.subset(["strength"+str(i) for i in range(12)])
    self.hdmg  = self.huse.subset(["whatever"])

    self.fixed = history.History()

  def test_hist(self):
    test = history.History()
    self.model.populate_history(test)
    self.assertEqual(test.size, 1)
    self.assertEqual(test.items, ["whatever"])
    self.model.init_history(test)
    self.assertAlmostEqual(test.get_scalar("whatever"), 0.0)

  def test_projection(self):
    P = self.model.projection(self.S, self.huse, self.Q, self.L, self.sliprule,
        self.T)
    self.assertEqual(P, tensors.SymSymR4.id())

  def test_damage_rate(self):
    rate = self.model.damage_rate(self.S, self.huse, self.Q, self.L,
        self.sliprule, self.T, self.fixed)
    self.assertEqual(rate.size, 1)
    self.assertEqual(rate.items, ["whatever"])
    self.assertAlmostEqual(rate.get_scalar("whatever"), 0)

class CommonSlipDamage:
  def test_damage_shear(self):
    act = self.model.d_damage_rate_d_shear(self.shear, self.sliprate,
        self.normal, self.damage)
    
    dfn = lambda s: self.model.damage_rate(list(s), self.sliprate, self.normal,
        self.damage)
    res = differentiate(dfn, np.array(self.shear))[0]

    self.assertTrue(np.allclose(act,res))

  def test_damage_slip(self):
    act = self.model.d_damage_rate_d_slip(self.shear, self.sliprate,
        self.normal, self.damage)
    
    dfn = lambda s: self.model.damage_rate(self.shear, list(s), self.normal,
        self.damage)
    res = differentiate(dfn, np.array(self.sliprate))[0]

    self.assertTrue(np.allclose(act,res))

  def test_damage_normal(self):
    act = self.model.d_damage_rate_d_normal(self.shear, self.sliprate,
        self.normal, self.damage)
    res = differentiate(lambda n: self.model.damage_rate(self.shear,
      self.sliprate, n, self.damage), self.normal)
    self.assertAlmostEqual(act, res)

  def test_damage_damage(self):
    act = self.model.d_damage_rate_d_damage(self.shear, self.sliprate,
        self.normal, self.damage)
    res = differentiate(lambda d: self.model.damage_rate(self.shear,
      self.sliprate, self.normal, d), self.damage)
    self.assertAlmostEqual(act, res)
    
class TestWorkSlipDamage(CommonSlipDamage, unittest.TestCase):
  def setUp(self):
    self.shear = [200.0,400.0,-400.0,150.0]
    self.sliprate = [1.0e-3,5.0e-2,-1.0e-1,1.0e-3]
    self.normal = 100.0
    self.damage = 0.5

    self.model = crystaldamage.WorkPlaneDamage()

  def test_damage_rate(self):
    drate = self.model.damage_rate(self.shear, self.sliprate, self.normal,
        self.damage)
    
    should = np.sum(np.array(self.shear) * np.array(self.sliprate))

    self.assertAlmostEqual(drate, should)

class CommonTransferFunctions:
  def test_map_damage(self):
    for v in self.vals:
      exact = self.function.d_map_d_damage(v, self.n)
      num = differentiate(lambda x: self.function.map(x, self.n), v)
      self.assertAlmostEqual(exact, num)

  def test_map_normal(self):
    for v in self.vals:
      exact = self.function.d_map_d_normal(v, self.n)
      num = differentiate(lambda x: self.function.map(v, x), self.n)
      self.assertAlmostEqual(exact, num)

class TestSigmoidTransfer(CommonTransferFunctions, unittest.TestCase):
  def setUp(self):
    self.c = 70.0
    self.beta = 2.1
    
    self.vals = np.array([1.0e-4,0.5*self.c,self.c,1.2*self.c])

    self.function = crystaldamage.SigmoidTransformation(self.c, self.beta)

    self.n = 100.0

  def test_value(self):
    x = np.array([self.function.map(v, self.n) for v in self.vals])
    y = np.piecewise(self.vals,
        [self.vals < self.c, self.vals >= self.c],
        [
          lambda x: 1.0/(1.0+(x/(self.c-x))**(-self.beta)),
          lambda x: 0.0*x+1.0
        ])

    self.assertTrue(np.allclose(x, y))

class TestSigmoidTransferP(CommonTransferFunctions, unittest.TestCase):
  def setUp(self):
    self.c = 70.0
    self.beta = 2.1
    
    self.vals = np.array([1.0e-4,0.5*self.c,self.c,1.2*self.c])

    self.function1 = crystaldamage.SigmoidTransformation(self.c, self.beta)
    self.function = crystaldamage.SwitchTransformation(self.function1)

    self.n = 100.0

  def test_value(self):
    x = np.array([self.function.map(v, self.n) for v in self.vals])
    y = np.piecewise(self.vals,
        [self.vals < self.c, self.vals >= self.c],
        [
          lambda x: 1.0/(1.0+(x/(self.c-x))**(-self.beta)),
          lambda x: 0.0*x+1.0
        ])

    self.assertTrue(np.allclose(x, y))

class TestSigmoidTransferN(CommonTransferFunctions, unittest.TestCase):
  def setUp(self):
    self.c = 70.0
    self.beta = 2.1
    
    self.vals = np.array([1.0e-4,0.5*self.c,self.c,1.2*self.c])

    self.function1 = crystaldamage.SigmoidTransformation(self.c, self.beta)
    self.function = crystaldamage.SwitchTransformation(self.function1)

    self.n = -1

  def test_value(self):
    x = np.array([self.function.map(v, self.n) for v in self.vals])
    y = np.array([0 for v in self.vals])

    self.assertTrue(np.allclose(x, y))
