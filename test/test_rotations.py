#!/usr/bin/env python

from neml.math import rotations

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

class ScalarMath(object):
  def test_scalar_division(self):
    self.assertTrue(np.allclose((self.A / self.s).quat, self.a / self.s))

  def test_scalar_multiply(self):
    self.assertTrue(np.allclose(self.a * self.s, (self.s * self.A).quat))
    self.assertTrue(np.allclose(self.a * self.s, (self.A * self.s).quat))

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

