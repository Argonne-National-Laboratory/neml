#!/usr/bin/env python

from neml.math import tensors

import unittest
import numpy as np
import numpy.linalg as la

class TestVector(unittest.TestCase):
  def setUp(self):
    self.a = np.array([2.2,-1.2,2.5])
    self.b = np.array([0.0,5.8,1.1])

    self.va = tensors.Vector(self.a)
    self.vb = tensors.Vector(self.b)

    self.s = 2.1
  
  def test_assign(self):
    self.vb = self.va
    self.assertTrue(np.allclose(self.vb.data, self.a))

  def test_norm(self):
    self.assertAlmostEqual(self.va.norm(), la.norm(self.a))

  def test_dot(self):
    self.assertAlmostEqual(np.dot(self.a,self.b), self.va.dot(self.vb))

  def test_smultiply(self):
    self.assertTrue(np.allclose(
      (self.s * self.va).data,
      self.s * self.a))
    self.assertTrue(np.allclose(
      (self.va * self.s).data,
      self.s * self.a))

  def test_sdivide(self):
    self.assertTrue(np.allclose(
      (self.va / self.s).data,
      self.a / self.s))

  def test_add(self):
    self.assertTrue(np.allclose(
      (self.va + self.vb).data,
      self.a + self.b))
    self.assertTrue(np.allclose(
      (self.vb + self.va).data,
      self.a + self.b))

  def test_subtract(self):
    self.assertTrue(np.allclose(
      (self.va - self.vb).data,
      self.a - self.b))
    self.assertTrue(np.allclose(
      (self.vb - self.va).data,
      self.b - self.a))

  def test_negate(self):
    self.assertTrue(np.allclose(
      (-self.va).data,
      -self.a))
    self.assertTrue(np.allclose(
      self.va.opposite().data,
      -self.a))

  def test_normalize(self):
    self.va.normalize()
    self.assertTrue(np.allclose(
      self.va.data, self.a / la.norm(self.a)))

  def test_cross(self):
    self.assertTrue(np.allclose(
      self.va.cross(self.vb).data,
      np.cross(self.a, self.b)))

  def test_equality(self):
    self.assertTrue(self.va == self.va)
    self.assertFalse(self.va == self.vb)

  def test_inequality(self):
    self.assertFalse(self.va != self.va)
    self.assertTrue(self.va != self.vb)

  def test_get(self):
    for i in range(3):
      self.assertTrue(np.isclose(self.va[i], self.a[i]))

  def test_set(self):
    for i in range(3):
      self.va[i] = 2.0
      self.a[i] = 2.0
      self.assertTrue(np.isclose(self.va[i], self.a[i]))

  def test_outer(self):
    self.assertEqual(self.va.outer(self.vb), tensors.outer(self.va, self.vb))
    self.assertEqual(tensors.RankTwo(np.outer(self.a, self.b)), 
        tensors.outer(self.va, self.vb))

class TestTensor(unittest.TestCase):
  def setUp(self):
    self.A = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.TA = tensors.RankTwo(self.A)
    self.B = np.array([[10.2,-9.3,2.5],[0.1,3.1,2.8],[0.1,3.2,-6.1]])
    self.TB = tensors.RankTwo(self.B)

    self.a = np.array([2.2,-1.2,2.5])
    self.va = tensors.Vector(self.a)

    self.s = 2.1

  def test_equality(self):
    self.assertEqual(self.TA, self.TA)

  def test_inequality(self):
    self.assertNotEqual(self.TA, self.TB)

  def test_get(self):
    self.assertTrue(np.isclose(self.TA[0,0], self.A[0,0]))

  def test_set(self):
    self.A[0,0] = 1.5
    self.assertTrue(np.isclose(self.A[0,0], 1.5))

  def test_scalar_mult(self):
    self.assertEqual(tensors.RankTwo(self.s*self.A), self.s * self.TA)
    self.assertEqual(tensors.RankTwo(self.A / self.s), self.TA / self.s)

  def test_add(self):
    self.assertEqual(tensors.RankTwo(self.A + self.B), self.TA + self.TB)
    self.assertEqual(tensors.RankTwo(self.A - self.B), self.TA - self.TB)

  def test_matrix_vector(self):
    self.assertEqual(tensors.Vector(np.dot(self.A, self.a)), self.TA*self.va)
    self.assertEqual(tensors.Vector(np.dot(self.a, self.A)), self.va*self.TA)

  def test_matrix_matrix(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.A, self.B)), self.TA*self.TB)

  def test_inverse(self):
    self.assertEqual(tensors.RankTwo(la.inv(self.A)), self.TA.inverse())

  def test_transpose(self):
    self.assertEqual(tensors.RankTwo(self.A.T), self.TA.transpose())

class TestSymmetric(unittest.TestCase):
  def setUp(self):
    self.A = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.A = 0.5*(self.A + self.A.T)
    self.TA = tensors.Symmetric(self.A)
    self.B = np.array([[10.2,-9.3,2.5],[0.1,3.1,2.8],[0.1,3.2,-6.1]])
    self.B = 0.5*(self.B + self.B.T)
    self.TB = tensors.Symmetric(self.B)

    self.a = np.array([2.2,-1.2,2.5])
    self.va = tensors.Vector(self.a)

    self.s = 2.1

  def test_equality(self):
    self.assertEqual(self.TA, self.TA)

  def test_inequality(self):
    self.assertNotEqual(self.TA, self.TB)

  def test_scalar_mult(self):
    self.assertEqual(tensors.Symmetric(self.s*self.A), self.s * self.TA)
    self.assertEqual(tensors.Symmetric(self.A / self.s), self.TA / self.s)

  def test_add(self):
    self.assertEqual(tensors.Symmetric(self.A + self.B), self.TA + self.TB)
    self.assertEqual(tensors.Symmetric(self.A - self.B), self.TA - self.TB)

  def test_matrix_vector(self):
    self.assertEqual(tensors.Vector(np.dot(self.A, self.a)), self.TA*self.va)
    self.assertEqual(tensors.Vector(np.dot(self.a, self.A)), self.va*self.TA)

  def test_matrix_matrix(self):
    self.assertEqual(tensors.Symmetric(np.dot(self.A, self.B)), self.TA*self.TB)

  def test_inverse(self):
    self.assertEqual(tensors.Symmetric(la.inv(self.A)), self.TA.inverse())

  def test_transpose(self):
    self.assertEqual(tensors.Symmetric(self.A.T), self.TA.transpose())

# Test various combinations of tensors
class TestComboTensorProducts(unittest.TestCase):
  def setUp(self):
    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)
    self.B = np.array([[10.2,-9.3,2.5],[0.1,3.1,2.8],[0.1,3.2,-6.1]])
    self.TB = tensors.RankTwo(self.B)

  def test_sym_general(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.S, self.B)), self.TS * self.TB)

  def test_general_sym(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.B, self.S)), self.TB * self.TS)
    
