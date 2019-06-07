#!/usr/bin/env python

from neml.math import tensors

import common

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

class TestSkew(unittest.TestCase):
  def setUp(self):
    self.A = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.A = 0.5*(self.A - self.A.T)
    self.TA = tensors.Skew(self.A)
    self.B = np.array([[10.2,-9.3,2.5],[0.1,3.1,2.8],[0.1,3.2,-6.1]])
    self.B = 0.5*(self.B - self.B.T)
    self.TB = tensors.Skew(self.B)

    self.a = np.array([2.2,-1.2,2.5])
    self.va = tensors.Vector(self.a)

    self.s = 2.1

  def test_equality(self):
    self.assertEqual(self.TA, self.TA)

  def test_inequality(self):
    self.assertNotEqual(self.TA, self.TB)

  def test_scalar_mult(self):
    self.assertEqual(tensors.Skew(self.s*self.A), self.s * self.TA)
    self.assertEqual(tensors.Skew(self.A / self.s), self.TA / self.s)

  def test_add(self):
    self.assertEqual(tensors.Skew(self.A + self.B), self.TA + self.TB)
    self.assertEqual(tensors.Skew(self.A - self.B), self.TA - self.TB)

  def test_matrix_vector(self):
    self.assertEqual(tensors.Vector(np.dot(self.A, self.a)), self.TA*self.va)
    self.assertEqual(tensors.Vector(np.dot(self.a, self.A)), self.va*self.TA)

  def test_matrix_matrix(self):
    self.assertEqual(tensors.Skew(np.dot(self.A, self.B)), self.TA*self.TB)

  def test_transpose(self):
    self.assertEqual(tensors.Skew(self.A.T), self.TA.transpose())

# Test various multiplicative combinations of tensors
class TestComboTensorMultiply(unittest.TestCase):
  def setUp(self):
    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)
    self.G = np.array([[10.2,-9.3,2.5],[0.1,3.1,2.8],[0.1,3.2,-6.1]])
    self.TG = tensors.RankTwo(self.G)
    self.W = np.array([[-5.0,7.1,1.0],[-0.2,0.25,1.2],[-0.4,0.4,-2]])
    self.W = 0.5*(self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

  def test_sym_general(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.S, self.G)), self.TS * self.TG)

  def test_general_sym(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.G, self.S)), self.TG * self.TS)

  def test_skew_general(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.W, self.G)), self.TW * self.TG)

  def test_general_skew(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.G, self.W)), self.TG * self.TW)

  def test_skew_sym(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.W, self.S)), self.TW * self.TS)
  
  def test_sym_skew(self):
    self.assertEqual(tensors.RankTwo(np.dot(self.S, self.W)), self.TS * self.TW)
   
class TestComboTensorAdd(unittest.TestCase):
  def setUp(self):
    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)
    self.G = np.array([[10.2,-9.3,2.5],[0.1,3.1,2.8],[0.1,3.2,-6.1]])
    self.TG = tensors.RankTwo(self.G)
    self.W = np.array([[-5.0,7.1,1.0],[-0.2,0.25,1.2],[-0.4,0.4,-2]])
    self.W = 0.5*(self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

  def test_add_sym_general(self):
    self.assertEqual(tensors.RankTwo(self.S + self.G), self.TS + self.TG)

  def test_add_general_sym(self):
    self.assertEqual(tensors.RankTwo(self.G + self.S), self.TG + self.TS)

  def test_add_skew_general(self):
    self.assertEqual(tensors.RankTwo(self.W + self.G), self.TW + self.TG)

  def test_add_general_skew(self):
    self.assertEqual(tensors.RankTwo(self.G + self.W), self.TG + self.TW)

  def test_sub_sym_general(self):
    self.assertEqual(tensors.RankTwo(self.S - self.G), self.TS - self.TG)

  def test_sub_general_sym(self):
    self.assertEqual(tensors.RankTwo(self.G - self.S), self.TG - self.TS)

  def test_sub_skew_general(self):
    self.assertEqual(tensors.RankTwo(self.W - self.G), self.TW - self.TG)

  def test_sub_general_skew(self):
    self.assertEqual(tensors.RankTwo(self.G - self.W), self.TG - self.TW)

class TestSymSym(unittest.TestCase):
  def setUp(self):
    self.SS1 = np.array([
      [ 5.99159801, -2.24342348,  0.26667281, -0.95466199,  3.98931478, -0.10846981],
      [ 1.86468226, -4.32391908, -7.82738638, -7.45008989,  5.89874777, 0.45820648],
      [-5.92565398,  2.4862829 , -6.02112389,  6.75455965,  4.65183463, 9.96900579],
      [ 0.60378883, -3.72189328, -7.63388446, -5.76559403, -0.3119789 , -1.1527258 ],
      [ 4.56813135, -6.06783828, -6.18341368,  8.06169686, -9.56928844, 9.08114655],
      [-8.25516614,  6.30663846,  7.2084381 , -7.38280703, -5.96279902, 8.9935982 ]])
    self.TSS1 = tensors.SymSym(self.SS1)

    self.SS2 = np.array([
      [-3.83767383, -8.63726504, -4.52095938,  9.35252323,  2.12800902, 3.26478511],
      [ 0.41705962,  3.95885105, -4.21676978,  4.12817198,  7.38839962, 5.79308578],
      [ 6.09635931,  2.31981366, -4.40237946, -5.51856189,  5.63572381, -5.55192385],
      [-0.97547288, -6.35708101, -4.35087656, -2.56567326,  4.32627031, 5.99408963],
      [ 6.30359707,  5.72926973,  2.47121354, -7.26333416, -5.08412215, -9.21872687],
      [-6.10780884,  1.01881487, -1.93491321,  6.13272186, -8.8721007, -2.97045116]])
    self.TSS2 = tensors.SymSym(self.SS2)

    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)

    self.scalar = 5.2

  def test_add(self):
    self.assertEqual(tensors.SymSym(self.SS1 + self.SS2), self.TSS2 + self.TSS1)
    self.assertEqual(tensors.SymSym(self.SS1 - self.SS2), self.TSS1 - self.TSS2)

  def test_equality(self):
    self.assertEqual(self.TSS1, self.TSS1)

  def test_inequality(self):
    self.assertNotEqual(self.TSS1, self.TSS2)

  def test_negate(self):
    self.assertEqual(tensors.SymSym(-self.SS1), -self.TSS1)

  def test_scalar_mult(self):
    self.assertEqual(tensors.SymSym(self.scalar * self.SS1), self.scalar * self.TSS1)
    self.assertEqual(tensors.SymSym(self.scalar * self.SS2), self.TSS2 * self.scalar)
    self.assertEqual(tensors.SymSym(self.SS1 / self.scalar), self.TSS1 / self.scalar)

  def test_product_sym_sym(self):
    self.assertEqual(tensors.SymSym(np.dot(self.SS1, self.SS2)), self.TSS1 * self.TSS2)

  def test_product_sym(self):
    self.assertEqual(tensors.Symmetric(common.usym(np.dot(self.SS1, common.sym(self.S)))), self.TSS1 * self.TS)


