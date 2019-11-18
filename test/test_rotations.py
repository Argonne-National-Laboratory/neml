#!/usr/bin/env python

from neml.math import rotations, tensors
from common import ms2ts

import unittest
import numpy as np
import numpy.linalg as la

def quat_mult(a, b):
  return np.concatenate((
    [a[0] * b[0] - np.dot(a[1:], b[1:])],
    a[0] * b[1:] + b[0] * a[1:] + np.cross(a[1:], b[1:])
    ))

def quat_exp(a):
  nv = la.norm(a[1:])

  nq = np.copy(a)
  a[0] = np.cos(nv)
  a[1:] *= np.sin(nv) / nv

  return a

class BasicMath(object):
  """
    Basic math operations
  """
  def test_quat(self):
    self.assertTrue(np.allclose(self.a, self.A.quat))

  def test_norm(self):
    self.assertAlmostEqual(la.norm(self.a), self.A.norm())

  def test_opposite(self):
    self.assertTrue(np.allclose(-self.a, -self.A.quat))

  def test_conjugate(self):
    aconj = np.copy(self.a)
    aconj[1:] = -aconj[1:]
    self.assertTrue(np.allclose(aconj, self.A.conj().quat))

  def test_flip(self):
    aflip = np.copy(self.a)
    aflip[0] = -aflip[0]

  def test_inverse(self):
    Ainv = self.A.inverse()
    iden = Ainv * self.A
    self.assertTrue(np.allclose(iden.quat, [1.0,0,0,0]))

  def test_quat_multiply(self):
    pv = quat_mult(self.a, self.b)
    PQ = self.A * self.B
    self.assertTrue(np.allclose(pv, PQ.quat))

  def test_quat_division(self):
    self.assertTrue(np.allclose(
      (self.A / self.B).quat,
      quat_mult(self.a, self.B.inverse().quat)))

  def test_exp(self):
    self.assertTrue(np.allclose((self.A.exp()).quat, quat_exp(self.a)))

  def test_dot(self):
    self.assertTrue(np.isclose(self.A.dot(self.B), 
      np.dot(self.A.quat, self.B.quat)))

  def test_matrix_composition(self):
    AB = self.A * self.B
    c = rotations.Quaternion(np.dot(self.A.to_product_matrix(), self.B.quat))
    self.assertTrue(np.allclose(AB.quat, c.quat))

class ScalarMath(object):
  def test_scalar_division(self):
    self.assertTrue(np.allclose((self.A / self.s).quat, self.a / self.s))

  def test_scalar_multiply(self):
    self.assertTrue(np.allclose(self.a * self.s, (self.s * self.A).quat))
    self.assertTrue(np.allclose(self.a * self.s, (self.A * self.s).quat))

class TestHash(unittest.TestCase):
  def setUp(self):
    self.Q1 = rotations.Quaternion(np.array([1.0,-2.1,3.21,4.0]))
    self.Q2 = rotations.Quaternion(np.array([5.0,-2.1,3.21,4.0]))
    self.Q3 = rotations.Quaternion(np.array([1.0000001,-2.1,3.21,4.0]))
    self.Q4 = rotations.Quaternion(np.array([4.0,-2.1,3.21,1.0]))

  def test_same(self):
    self.assertEqual(self.Q1.hash, self.Q1.hash)

  def test_not_same(self):
    self.assertNotEqual(self.Q1.hash, self.Q2.hash)
    self.assertNotEqual(self.Q1.hash, self.Q3.hash)
    self.assertNotEqual(self.Q1.hash, self.Q4.hash)

class TestBasicMathQuaternion(unittest.TestCase, BasicMath, ScalarMath):
  def setUp(self):
    self.a = np.array([2.0,-4.0,5.0,1.0])
    self.A = rotations.Quaternion(self.a)
    self.b = np.array([-1.0,3.0,2.0,0.5])
    self.B = rotations.Quaternion(self.b)

    self.s = 2.2

class TestBasicMathCOrientation(unittest.TestCase, BasicMath):
  def setUp(self):
    self.A = rotations.Orientation(31.0, 59.0, 80.0, angle_type = "degrees",
        convention = "bunge")
    self.a = np.copy(self.A.quat)
    n = np.array([1,2.0,3])
    n /= la.norm(n)
    self.B = rotations.Orientation(n, np.pi/3, angle_type = "radians")
    self.b = np.copy(self.B.quat)

class Conversion(object):
  def test_matrix(self):
    M = self.q.to_matrix()
    q2 = rotations.Orientation(M)
    self.assertTrue(np.allclose(self.q.quat, q2.quat) or 
        np.allclose(self.q.quat, -q2.quat))

  def test_axis_angle(self):
    (axis, angle) = self.q.to_axis_angle()
    q2 = rotations.Orientation(axis, angle)
    self.assertTrue(np.allclose(self.q.quat, q2.quat) or 
        np.allclose(self.q.quat, -q2.quat))

  def test_rodrigues(self):
    rv = self.q.to_rodrigues()
    q2 = rotations.Orientation(rv)
    self.assertTrue(np.allclose(self.q.quat, q2.quat) or 
        np.allclose(self.q.quat, -q2.quat))

  def test_euler_kocks(self):
    a, b, c = self.q.to_euler()
    q2 = rotations.Orientation(a, b, c)
    self.assertTrue(np.allclose(self.q.quat, q2.quat) or 
        np.allclose(self.q.quat, -q2.quat))

  def test_euler_bunge(self):
    a, b, c = self.q.to_euler(angle_type = "degrees", convention = "bunge")
    q2 = rotations.Orientation(a, b, c, angle_type = "degrees", 
        convention = "bunge")
    self.assertTrue(np.allclose(self.q.quat, q2.quat) or 
        np.allclose(self.q.quat, -q2.quat))

  def test_hopf(self):
    a, b, c = self.q.to_hopf(angle_type = "degrees")
    q2 = rotations.Orientation(a, b, c, angle_type = "degrees",
        convention = "hopf")
    self.assertTrue(np.allclose(self.q.quat, q2.quat) or 
        np.allclose(self.q.quat, -q2.quat))

  def test_hyperspherical(self):
    a1, a2, a3 = self.q.to_hyperspherical(angle_type = "degrees")
    q2 = rotations.Orientation(a1, a2, a3, angle_type = "degrees",
        convention = "hyperspherical")
    self.assertTrue(np.allclose(self.q.quat, q2.quat) or 
        np.allclose(self.q.quat, -q2.quat))

class TestEasyCaseConversion(unittest.TestCase, Conversion):
  def setUp(self):
    self.q = rotations.Orientation(30.0, 60.0, 80.0, angle_type = "degrees")

class TestSpecialCaseConversion(unittest.TestCase, Conversion):
  def setUp(self):
    self.q = rotations.Orientation(0.0, 0.0, 10.0, angle_type = "degrees")

class TestApplyThings(unittest.TestCase):
  def setUp(self):
    self.q = rotations.Orientation(30.0, 60.0, 80.0, angle_type = "degrees")
    self.Q = self.q.to_matrix()
    self.v = np.array([1.2,-2.0,3.0])
    self.Tv = tensors.Vector(self.v)

    self.A = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.TA = tensors.RankTwo(self.A)

    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)

    self.W = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.W = 0.5*(self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

    self.R1 = np.array([[[[ 7.09627147,  9.22330744, -1.36602973],
             [-7.86118175, -1.6342633 , -5.75516189],
             [ 2.61734248,  6.40678382,  3.37981603]],
            [[ 5.65100254, -7.88797059,  7.31396665],
             [-6.35471595,  5.67698069, -8.18795178],
             [ 9.10447016,  8.91183436, -6.65254333]],
            [[ 3.20429862,  2.99308849,  4.0035241 ],
             [-4.02440197, -4.39975872, -4.33542791],
             [ 9.36746226, -2.91156335,  4.51572032]]],
           [[[-9.23675199,  8.63546962,  6.83448027],
             [ 4.35044123,  2.24508666,  9.80054664],
             [ 0.30835223, -4.05208575,  5.68966326]],
            [[ 6.40300092, -8.25998136,  5.63566553],
             [-5.02801101,  5.64005224, -7.39586166],
             [ 5.90893633,  6.02074669,  1.37112738]],
            [[-2.68485216, -4.67660156,  3.52618441],
             [-2.52484812, -0.08561168,  3.39072868],
             [ 9.11295675,  2.63102786, -4.82285415]]],
           [[[ 8.31973154,  4.76081593,  4.38377207],
             [ 6.22896742, -3.83995097,  5.37501029],
             [-0.16770967,  7.9453854 , -4.95548491]],
            [[-5.67884611, -8.44970885, -7.42037867],
             [-5.19908193, -7.87006493,  1.65949787],
             [-3.25934672,  6.27340198,  5.98643056]],
            [[-4.20166968, -2.38276224,  3.04551936],
             [ 3.68445989, -5.84357996,  3.61183543],
             [ 1.54886677,  3.3659842 ,  6.43067337]]]])
    self.TR1 = tensors.RankFour(self.R1)

    self.SS1 = np.array([
      [ 5.99159801, -2.24342348,  0.26667281, -0.95466199,  3.98931478, -0.10846981],
      [ 1.86468226, -4.32391908, -7.82738638, -7.45008989,  5.89874777, 0.45820648],
      [-5.92565398,  2.4862829 , -6.02112389,  6.75455965,  4.65183463, 9.96900579],
      [ 0.60378883, -3.72189328, -7.63388446, -5.76559403, -0.3119789 , -1.1527258 ],
      [ 4.56813135, -6.06783828, -6.18341368,  8.06169686, -9.56928844, 9.08114655],
      [-8.25516614,  6.30663846,  7.2084381 , -7.38280703, -5.96279902, 8.9935982 ]])
    self.SSS1 = ms2ts(self.SS1)
    self.TSS1 = tensors.SymSymR4(self.SS1)

  def test_apply_vector(self):
    rv = self.q.apply(self.Tv)
    self.assertAlmostEqual(rv.norm(), self.Tv.norm())
    rrv = self.q.inverse().apply(rv)
    self.assertTrue(np.allclose(rrv.data, self.v.data))

    self.assertEqual(self.q.apply(self.Tv), tensors.Vector(np.dot(self.Q, self.v)))
  
  def test_apply_ranktwo(self):
    self.assertEqual(self.q.apply(self.TA), 
        tensors.RankTwo(np.dot(self.Q, np.dot(self.A, self.Q.T))))

  def test_apply_symmetric(self):
    self.assertEqual(self.q.apply(self.TS),
        tensors.Symmetric(np.dot(self.Q, np.dot(self.S, self.Q.T))))

  def test_apply_skew(self):
    self.assertEqual(self.q.apply(self.TW),
        tensors.Skew(np.dot(self.Q, np.dot(self.W, self.Q.T))))

  @staticmethod
  def rotate_fourth(T, Q):
    return np.einsum('ip,jq,pqrs,kr,ls', Q, Q, T, Q, Q)

  def test_apply_rankfour(self):
    self.assertEqual(self.q.apply(self.TR1), 
        tensors.RankFour(TestApplyThings.rotate_fourth(self.R1, self.Q)))

  def test_apply_symsym(self):
    self.assertEqual(self.q.apply(self.TSS1),
        tensors.RankFour(TestApplyThings.rotate_fourth(self.SSS1, self.Q)
          ).to_sym())

class TestExpIntegration(unittest.TestCase):
  def setUp(self):
    self.q = rotations.Orientation(30.0, 60.0, 80.0, angle_type = "degrees")

    self.W = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.W = 0.5*(self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

  def test_w2q(self):
    self.assertEqual(self.TW, rotations.wlog(rotations.wexp(self.TW)))

  def test_q2w(self):
    self.assertTrue(np.allclose(self.q.quat,
      rotations.wexp(rotations.wlog(self.q)).quat))

  def test_correctly_integrates(self):
    # So we are rotating about x where t is the rotation in radians
    R = lambda t: np.array([[1.0,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t), np.cos(t)]])
    q0 = rotations.Orientation(R(0.0))

    spin = tensors.Skew(np.array([
      [0,0,0],
      [0,0,-1.0],
      [0,1.0,0]]))
    
    nsteps = 1000 # Well it's not the most accurate integration
    ang = np.pi/4.0

    qf = rotations.Orientation(R(ang))
    
    q = q0
    for i in range(nsteps):
      q = rotations.wexp(spin * ang/nsteps) * q

    self.assertTrue(np.allclose(q.quat, qf.quat))

class TestFromVectors(unittest.TestCase):
  def setUp(self):
    self.x = tensors.Vector([0.0,0.0,1.0])
    self.y = tensors.Vector([1.0,0.0,0.0])
    self.q = rotations.Orientation(self.x,self.y)

  def test_sends_correctly(self):
    x0 = tensors.Vector([1.0,0,0])
    y0 = tensors.Vector([0.0,1.0,0.0])

    xp = self.q.apply(x0)
    self.assertEqual(self.x, xp)

    yp = self.q.apply(y0)
    self.assertEqual(self.y, yp)

class TestDistance(unittest.TestCase):
  def setUp(self):
    self.q = rotations.Orientation(30.0, 60.0, 80.0, angle_type = "degrees")

  def test_same_zero(self):
    self.assertTrue(np.isclose(self.q.distance(self.q), 0))
    self.assertTrue(np.isclose(rotations.distance(self.q,self.q), 0))

class TestOrientVectors(unittest.TestCase):
  def setUp(self):
    self.v1 = tensors.Vector([1.0,4,2])
    self.v1.normalize()
    self.v2 = tensors.Vector([2,1,-6])
    self.v2.normalize()

  def test_single_orient(self):
    q = rotations.rotate_to(self.v1, self.v2)
    v3 = q.apply(self.v1)
    self.assertEqual(self.v2, v3)

  def test_family_orient(self):
    for a in np.linspace(0,2*np.pi):
      q = rotations.rotate_to_family(self.v1, self.v2, a)
      v3 = q.apply(self.v1)
      self.assertEqual(self.v2,v3)

