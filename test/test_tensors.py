#!/usr/bin/env python

from neml.math import tensors

import common

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

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

class TestRankTwo(unittest.TestCase):
  def setUp(self):
    self.A = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.TA = tensors.RankTwo(self.A)
    self.B = np.array([[10.2,-9.3,2.5],[0.1,3.1,2.8],[0.1,3.2,-6.1]])
    self.TB = tensors.RankTwo(self.B)

    self.a = np.array([2.2,-1.2,2.5])
    self.va = tensors.Vector(self.a)

    self.s = 2.1

  def test_id(self):
    self.assertTrue(tensors.RankTwo.id(), tensors.RankTwo(np.eye(3)))

  def test_norm(self):
    self.assertTrue(np.isclose(self.TA.norm(), np.sqrt(np.sum(self.A*self.A))))

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

  def test_norm(self):
    self.assertTrue(np.isclose(self.TA.norm(), np.sqrt(np.sum(self.A*self.A))))

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

  def test_id(self):
    self.assertEqual(tensors.Symmetric.id(), tensors.Symmetric(np.eye(3)))

  def test_trace(self):
    self.assertTrue(np.isclose(self.TA.trace(), np.trace(self.A)))

  def test_dev(self):
    self.assertTrue(self.TA.dev(), tensors.Symmetric(
      self.A - np.trace(self.A)/3.0 * np.eye(3)))

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

  def test_contract_general_sym(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.G, self.S),
      self.TG.contract(self.TS)))

  def test_contract_sym_general(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.G, self.S),
      self.TS.contract(self.TG)))

  def test_contract_general_skew(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.G, self.W),
      self.TG.contract(self.TW)))

  def test_contract_skew_general(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.G, self.W),
      self.TW.contract(self.TG)))

  def test_contract_skew_sym(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.W, self.S),
      self.TW.contract(self.TS)))

  def test_contract_sym_skew(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.W, self.S),
      self.TS.contract(self.TW)))

  def test_contract_general_general(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.G, self.G),
      self.TG.contract(self.TG)))

  def test_contract_sym_sym(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.S, self.S),
      self.TS.contract(self.TS)))

  def test_contract_skew_skew(self):
    self.assertTrue(np.isclose(
      np.einsum('ij,ij', self.W, self.W),
      self.TW.contract(self.TW)))
   
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

class TestRankFour(unittest.TestCase):
  def setUp(self):
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

    self.R2 = np.array([[[[-8.03675620e+00,  2.58575052e+00,  2.44069661e+00],
             [ 4.75021663e+00,  1.24463394e+00, -8.69751301e-01],
             [-1.46310894e+00, -1.15053235e+00, -3.75342982e+00]],
            [[-7.64033956e+00,  4.19956720e+00, -4.87644982e+00],
             [ 1.06577507e+00,  8.94272637e+00,  6.57264250e-01],
             [-4.22613258e+00, -5.08830314e+00,  1.57718186e+00]],
            [[-4.02243082e+00, -4.75463781e+00, -8.88662152e+00],
             [-1.30383950e+00, -1.98063574e+00, -3.18963544e+00],
             [-7.52071674e+00,  1.08931933e+00,  2.86988431e+00]]],
           [[[ 5.28621060e+00, -6.83799668e+00,  8.98005935e+00],
             [-7.92741122e+00,  5.75699425e-01,  1.66782544e+00],
             [ 2.60041984e+00, -1.04476986e-02, -6.12424787e+00]],
            [[-3.73727368e+00,  6.59764771e+00, -1.18045587e+00],
             [ 4.08567441e+00,  2.66148943e+00, -6.82495588e-01],
             [-1.64417262e+00,  5.33119298e+00,  8.11045988e-03]],
            [[-5.90193883e+00, -2.63316107e+00,  5.61381825e+00],
             [-6.08591194e+00,  8.77285539e+00, -7.15230533e+00],
             [ 3.15093096e+00,  1.41350149e+00,  1.11702016e+00]]],
           [[[-9.61472764e-01, -1.91492497e+00,  9.48275324e+00],
             [ 6.68841134e+00,  3.23412041e+00, -3.41944541e+00],
             [-9.80203467e+00,  6.58425335e+00, -2.16548636e+00]],
            [[ 6.63950740e+00,  3.91551441e+00, -8.98229111e+00],
             [ 9.84606756e+00, -8.16145090e+00,  8.41929062e-01],
             [-1.93839620e+00,  7.44485127e+00, -2.70832414e+00]],
            [[ 9.79265531e+00, -1.18212395e+00, -5.39433704e+00],
             [ 4.87152614e+00,  9.47287450e+00,  5.53838514e+00],
             [ 9.30443367e+00,  1.27090319e+00,  1.60409739e+00]]]])
    self.TR2 = tensors.RankFour(self.R2)

    self.SS = np.array([
      [ 5.99159801, -2.24342348,  0.26667281, -0.95466199,  3.98931478, -0.10846981],
      [ 1.86468226, -4.32391908, -7.82738638, -7.45008989,  5.89874777, 0.45820648],
      [-5.92565398,  2.4862829 , -6.02112389,  6.75455965,  4.65183463, 9.96900579],
      [ 0.60378883, -3.72189328, -7.63388446, -5.76559403, -0.3119789 , -1.1527258 ],
      [ 4.56813135, -6.06783828, -6.18341368,  8.06169686, -9.56928844, 9.08114655],
      [-8.25516614,  6.30663846,  7.2084381 , -7.38280703, -5.96279902, 8.9935982 ]])
    self.SS_full = common.ms2ts(self.SS)
    self.TSS = tensors.SymSymR4(self.SS)

    self.SW = np.array([
      [ 5.43434005, -6.55983214,  0.29737664],
      [-4.77472172, -8.51287287, -3.19380185],
      [ 4.43407952, -6.02555614,  5.87786914],
      [ 1.89488869, -5.65383917,  8.83717547],
      [-7.18030867,  1.56100537, -9.83238641],
      [-4.52369317, -3.07284914, -7.54966999]])
    self.SW_full = common.ws2ts(self.SW)
    self.TSW = tensors.SymSkewR4(self.SW)

    self.WS = np.array([
      [-8.3567359 , -5.39728818, -8.00844442, -8.33365112, -0.97903364, -8.23943149],
      [-6.97125417,  4.34802055,  7.06281056, -1.57511617,  7.83359933, -9.37625432],
      [-6.0799489 , -6.0309543 ,  3.68575895,  8.84296976,  6.55799427, -9.22029379]])
    self.WS_full = common.wws2ts(self.WS)
    self.TWS = tensors.SkewSymR4(self.WS)

    self.scalar = -2.2

    self.G = np.array([[ 9.50640677,  1.79084726, -2.8877036 ],
       [-1.63159958,  2.52866904, -8.71585042],
       [ 5.01859685, -8.7324075 , -0.42919134]])
    self.TG = tensors.RankTwo(self.G)

    self.S = np.array([[ 6.19999242, -6.95811611, -6.02901899],
           [ 8.38508084,  6.01607694,  6.79839425],
           [-4.4214246 , -2.36795313, -8.84070728]])
    self.S = 0.5*(self.S+self.S.T)
    self.TS = tensors.Symmetric(self.S)

    self.W = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.W = 0.5 * (self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

  def test_add(self):
    self.assertEqual(tensors.RankFour(self.R1 + self.R2), self.TR2 + self.TR1)
    self.assertEqual(tensors.RankFour(self.R1 - self.R2), self.TR1 - self.TR2)

  def test_equality(self):
    self.assertEqual(self.TR1, self.TR1)

  def test_inequality(self):
    self.assertNotEqual(self.TR1, self.TR2)

  def test_negate(self):
    self.assertEqual(tensors.RankFour(-self.R1), -self.TR1)

  def test_scalar_mult(self):
    self.assertEqual(tensors.RankFour(self.scalar * self.R1), self.scalar * self.TR1)
    self.assertEqual(tensors.RankFour(self.scalar * self.R2), self.TR2 * self.scalar)
    self.assertEqual(tensors.RankFour(self.R1 / self.scalar), self.TR1 / self.scalar)
  
  def test_double_contraction(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.R1, self.R2)), self.TR1 * self.TR2)

  def test_sym_sym(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.R1, self.SS_full)), self.TR1 * self.TSS)

  def test_sym_skew(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.R1, self.SW_full)), self.TR1 * self.TSW)

  def test_skew_sym(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.R1, self.WS_full)), self.TR1 * self.TWS)

  def test_ranktwo(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.R1, self.G)), self.TR1 * self.TG)

  def test_symmetric(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.R1, self.S)), self.TR1 * self.TS)

  def test_skew(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.R1, self.W)), self.TR1 * self.TW)

  def test_get(self):
    self.assertTrue(np.isclose(self.R1[1,2,0,1], self.TR1[1,2,0,1]))

  def test_set(self):
    self.TR1[1,1,1,1] = 4.0
    self.assertTrue(np.isclose(self.TR1[1,1,1,1],4.0))

class TestSymSymR4(unittest.TestCase):
  def setUp(self):
    self.SS1 = np.array([
      [ 5.99159801, -2.24342348,  0.26667281, -0.95466199,  3.98931478, -0.10846981],
      [ 1.86468226, -4.32391908, -7.82738638, -7.45008989,  5.89874777, 0.45820648],
      [-5.92565398,  2.4862829 , -6.02112389,  6.75455965,  4.65183463, 9.96900579],
      [ 0.60378883, -3.72189328, -7.63388446, -5.76559403, -0.3119789 , -1.1527258 ],
      [ 4.56813135, -6.06783828, -6.18341368,  8.06169686, -9.56928844, 9.08114655],
      [-8.25516614,  6.30663846,  7.2084381 , -7.38280703, -5.96279902, 8.9935982 ]])
    self.SS1_full = common.ms2ts(self.SS1)
    self.TSS1 = tensors.SymSymR4(self.SS1)

    self.SS2 = np.array([
      [-3.83767383, -8.63726504, -4.52095938,  9.35252323,  2.12800902, 3.26478511],
      [ 0.41705962,  3.95885105, -4.21676978,  4.12817198,  7.38839962, 5.79308578],
      [ 6.09635931,  2.31981366, -4.40237946, -5.51856189,  5.63572381, -5.55192385],
      [-0.97547288, -6.35708101, -4.35087656, -2.56567326,  4.32627031, 5.99408963],
      [ 6.30359707,  5.72926973,  2.47121354, -7.26333416, -5.08412215, -9.21872687],
      [-6.10780884,  1.01881487, -1.93491321,  6.13272186, -8.8721007, -2.97045116]])
    self.TSS2 = tensors.SymSymR4(self.SS2)

    self.SW = np.array([
      [ 5.43434005, -6.55983214,  0.29737664],
      [-4.77472172, -8.51287287, -3.19380185],
      [ 4.43407952, -6.02555614,  5.87786914],
      [ 1.89488869, -5.65383917,  8.83717547],
      [-7.18030867,  1.56100537, -9.83238641],
      [-4.52369317, -3.07284914, -7.54966999]])
    self.SW_full = common.ws2ts(self.SW)
    self.TSW = tensors.SymSkewR4(self.SW)

    self.WS = np.array([
      [-8.3567359 , -5.39728818, -8.00844442, -8.33365112, -0.97903364, -8.23943149],
      [-6.97125417,  4.34802055,  7.06281056, -1.57511617,  7.83359933, -9.37625432],
      [-6.0799489 , -6.0309543 ,  3.68575895,  8.84296976,  6.55799427, -9.22029379]])
    self.WS_full = common.wws2ts(self.WS)
    self.TWS = tensors.SkewSymR4(self.WS)

    self.R = np.array([[[[-8.03675620e+00,  2.58575052e+00,  2.44069661e+00],
             [ 4.75021663e+00,  1.24463394e+00, -8.69751301e-01],
             [-1.46310894e+00, -1.15053235e+00, -3.75342982e+00]],
            [[-7.64033956e+00,  4.19956720e+00, -4.87644982e+00],
             [ 1.06577507e+00,  8.94272637e+00,  6.57264250e-01],
             [-4.22613258e+00, -5.08830314e+00,  1.57718186e+00]],
            [[-4.02243082e+00, -4.75463781e+00, -8.88662152e+00],
             [-1.30383950e+00, -1.98063574e+00, -3.18963544e+00],
             [-7.52071674e+00,  1.08931933e+00,  2.86988431e+00]]],
           [[[ 5.28621060e+00, -6.83799668e+00,  8.98005935e+00],
             [-7.92741122e+00,  5.75699425e-01,  1.66782544e+00],
             [ 2.60041984e+00, -1.04476986e-02, -6.12424787e+00]],
            [[-3.73727368e+00,  6.59764771e+00, -1.18045587e+00],
             [ 4.08567441e+00,  2.66148943e+00, -6.82495588e-01],
             [-1.64417262e+00,  5.33119298e+00,  8.11045988e-03]],
            [[-5.90193883e+00, -2.63316107e+00,  5.61381825e+00],
             [-6.08591194e+00,  8.77285539e+00, -7.15230533e+00],
             [ 3.15093096e+00,  1.41350149e+00,  1.11702016e+00]]],
           [[[-9.61472764e-01, -1.91492497e+00,  9.48275324e+00],
             [ 6.68841134e+00,  3.23412041e+00, -3.41944541e+00],
             [-9.80203467e+00,  6.58425335e+00, -2.16548636e+00]],
            [[ 6.63950740e+00,  3.91551441e+00, -8.98229111e+00],
             [ 9.84606756e+00, -8.16145090e+00,  8.41929062e-01],
             [-1.93839620e+00,  7.44485127e+00, -2.70832414e+00]],
            [[ 9.79265531e+00, -1.18212395e+00, -5.39433704e+00],
             [ 4.87152614e+00,  9.47287450e+00,  5.53838514e+00],
             [ 9.30443367e+00,  1.27090319e+00,  1.60409739e+00]]]])
    self.TR = tensors.RankFour(self.R)

    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)

    self.S2 = np.array([[10.2,-9.3,2.5],[0.1,3.1,2.8],[0.1,3.2,-6.1]])
    self.S2 = 0.5*(self.S2 + self.S2.T)
    self.TS2 = tensors.Symmetric(self.S2)

    self.scalar = 5.2

    self.G = np.array([[ 9.50640677,  1.79084726, -2.8877036 ],
       [-1.63159958,  2.52866904, -8.71585042],
       [ 5.01859685, -8.7324075 , -0.42919134]])
    self.TG = tensors.RankTwo(self.G)

    self.W = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.W = 0.5 * (self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

  def test_to_full(self):
    full_np = common.ms2ts(self.SS1)
    full_t = tensors.RankFour(full_np)
    full = self.TSS1.to_full()
    self.assertEqual(full_t, full)

  def test_from_full(self):
    full = self.TSS1.to_full()
    new = full.to_sym()
    self.assertEqual(self.TSS1, new)

  def test_add(self):
    self.assertEqual(tensors.SymSymR4(self.SS1 + self.SS2), self.TSS2 + self.TSS1)
    self.assertEqual(tensors.SymSymR4(self.SS1 - self.SS2), self.TSS1 - self.TSS2)

  def test_equality(self):
    self.assertEqual(self.TSS1, self.TSS1)

  def test_inequality(self):
    self.assertNotEqual(self.TSS1, self.TSS2)

  def test_negate(self):
    self.assertEqual(tensors.SymSymR4(-self.SS1), -self.TSS1)

  def test_scalar_mult(self):
    self.assertEqual(tensors.SymSymR4(self.scalar * self.SS1), self.scalar * self.TSS1)
    self.assertEqual(tensors.SymSymR4(self.scalar * self.SS2), self.TSS2 * self.scalar)
    self.assertEqual(tensors.SymSymR4(self.SS1 / self.scalar), self.TSS1 / self.scalar)

  def test_product_sym_sym(self):
    self.assertEqual(tensors.SymSymR4(np.dot(self.SS1, self.SS2)), self.TSS1 * self.TSS2)

  def test_product_sym_full(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.SS1_full, self.R)), self.TSS1 * self.TR)

  def test_product_sym_skew(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.SS1_full, self.SW_full)), self.TSS1 * self.TSW)

  def test_product_skew_sym(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.SS1_full, self.WS_full)), self.TSS1 * self.TWS)

  def test_product_symmetric(self):
    self.assertEqual(tensors.Symmetric(common.usym(np.dot(self.SS1, common.sym(self.S)))), self.TSS1 * self.TS)

  def test_product_general(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.SS1_full, self.G)), self.TSS1 * self.TG)

  def test_product_skew(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.SS1_full, self.W)), self.TSS1 * self.TW)

  def test_douter(self):
    self.assertEqual(tensors.SymSymR4(common.ts2ms(np.einsum('ij,kl', self.S, self.S2))), tensors.douter(self.TS, self.TS2))

  def test_id(self):
    id_t = tensors.SymSymR4(np.eye(6))
    self.assertEqual(id_t, tensors.SymSymR4.id())

  def test_id_dev(self):
    ot = np.zeros((6,6))
    ot[:3,:3] = 1.0/3.0
    id_t = tensors.SymSymR4(np.eye(6) - ot)
    self.assertEqual(id_t, tensors.SymSymR4.id_dev())

class TestSymSkewR4(unittest.TestCase):
  def setUp(self):
    self.SW1 = np.array([
      [ 5.43434005, -6.55983214,  0.29737664],
      [-4.77472172, -8.51287287, -3.19380185],
      [ 4.43407952, -6.02555614,  5.87786914],
      [ 1.89488869, -5.65383917,  8.83717547],
      [-7.18030867,  1.56100537, -9.83238641],
      [-4.52369317, -3.07284914, -7.54966999]])
    self.SW1_full = common.ws2ts(self.SW1)
    self.TSW1 = tensors.SymSkewR4(self.SW1)

    self.SW2 = np.array([
      [ 7.90885123, -1.89089468, -6.95528566],
      [-2.53495619,  9.47533071, -2.76302205],
      [-8.57887706,  4.21216331, -7.68619983],
      [-5.45955495,  2.0523769 , -9.71153458],
      [-5.61696943, -4.02142773, -6.41654212],
      [-8.76272792, -3.60354692,  2.7402794 ]])

    self.SW2_full = common.ws2ts(self.SW2)
    self.TSW2 = tensors.SymSkewR4(self.SW2)

    self.WS = np.array([
      [-8.3567359 , -5.39728818, -8.00844442, -8.33365112, -0.97903364, -8.23943149],
      [-6.97125417,  4.34802055,  7.06281056, -1.57511617,  7.83359933, -9.37625432],
      [-6.0799489 , -6.0309543 ,  3.68575895,  8.84296976,  6.55799427, -9.22029379]])
    self.WS_full = common.wws2ts(self.WS)
    self.TWS = tensors.SkewSymR4(self.WS)

    self.SS = np.array([
      [ 5.99159801, -2.24342348,  0.26667281, -0.95466199,  3.98931478, -0.10846981],
      [ 1.86468226, -4.32391908, -7.82738638, -7.45008989,  5.89874777, 0.45820648],
      [-5.92565398,  2.4862829 , -6.02112389,  6.75455965,  4.65183463, 9.96900579],
      [ 0.60378883, -3.72189328, -7.63388446, -5.76559403, -0.3119789 , -1.1527258 ],
      [ 4.56813135, -6.06783828, -6.18341368,  8.06169686, -9.56928844, 9.08114655],
      [-8.25516614,  6.30663846,  7.2084381 , -7.38280703, -5.96279902, 8.9935982 ]])
    self.SS_full = common.ms2ts(self.SS)
    self.TSS = tensors.SymSymR4(self.SS)

    self.R = np.array([[[[-8.03675620e+00,  2.58575052e+00,  2.44069661e+00],
             [ 4.75021663e+00,  1.24463394e+00, -8.69751301e-01],
             [-1.46310894e+00, -1.15053235e+00, -3.75342982e+00]],
            [[-7.64033956e+00,  4.19956720e+00, -4.87644982e+00],
             [ 1.06577507e+00,  8.94272637e+00,  6.57264250e-01],
             [-4.22613258e+00, -5.08830314e+00,  1.57718186e+00]],
            [[-4.02243082e+00, -4.75463781e+00, -8.88662152e+00],
             [-1.30383950e+00, -1.98063574e+00, -3.18963544e+00],
             [-7.52071674e+00,  1.08931933e+00,  2.86988431e+00]]],
           [[[ 5.28621060e+00, -6.83799668e+00,  8.98005935e+00],
             [-7.92741122e+00,  5.75699425e-01,  1.66782544e+00],
             [ 2.60041984e+00, -1.04476986e-02, -6.12424787e+00]],
            [[-3.73727368e+00,  6.59764771e+00, -1.18045587e+00],
             [ 4.08567441e+00,  2.66148943e+00, -6.82495588e-01],
             [-1.64417262e+00,  5.33119298e+00,  8.11045988e-03]],
            [[-5.90193883e+00, -2.63316107e+00,  5.61381825e+00],
             [-6.08591194e+00,  8.77285539e+00, -7.15230533e+00],
             [ 3.15093096e+00,  1.41350149e+00,  1.11702016e+00]]],
           [[[-9.61472764e-01, -1.91492497e+00,  9.48275324e+00],
             [ 6.68841134e+00,  3.23412041e+00, -3.41944541e+00],
             [-9.80203467e+00,  6.58425335e+00, -2.16548636e+00]],
            [[ 6.63950740e+00,  3.91551441e+00, -8.98229111e+00],
             [ 9.84606756e+00, -8.16145090e+00,  8.41929062e-01],
             [-1.93839620e+00,  7.44485127e+00, -2.70832414e+00]],
            [[ 9.79265531e+00, -1.18212395e+00, -5.39433704e+00],
             [ 4.87152614e+00,  9.47287450e+00,  5.53838514e+00],
             [ 9.30443367e+00,  1.27090319e+00,  1.60409739e+00]]]])
    self.TR = tensors.RankFour(self.R)

    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)

    self.scalar = 5.2

    self.G = np.array([[ 9.50640677,  1.79084726, -2.8877036 ],
       [-1.63159958,  2.52866904, -8.71585042],
       [ 5.01859685, -8.7324075 , -0.42919134]])
    self.TG = tensors.RankTwo(self.G)

    self.W = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.W = 0.5 * (self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

  def test_to_full(self):
    full_np = common.ws2ts(self.SW1)
    full_t = tensors.RankFour(full_np)
    full = self.TSW1.to_full()
    self.assertEqual(full_t, full)

  def test_from_full(self):
    full = self.TSW1.to_full()
    new = full.to_symskew()
    self.assertEqual(self.TSW1, new)

  def test_add(self):
    self.assertEqual(tensors.SymSkewR4(self.SW1 + self.SW2), self.TSW2 + self.TSW1)
    self.assertEqual(tensors.SymSkewR4(self.SW1 - self.SW2), self.TSW1 - self.TSW2)

  def test_equality(self):
    self.assertEqual(self.TSW1, self.TSW1)

  def test_inequality(self):
    self.assertNotEqual(self.TSW1, self.TSW2)

  def test_negate(self):
    self.assertEqual(tensors.SymSkewR4(-self.SW1), -self.TSW1)

  def test_scalar_mult(self):
    self.assertEqual(tensors.SymSkewR4(self.scalar * self.SW1), self.scalar * self.TSW1)
    self.assertEqual(tensors.SymSkewR4(self.scalar * self.SW2), self.TSW2 * self.scalar)
    self.assertEqual(tensors.SymSkewR4(self.SW1 / self.scalar), self.TSW1 / self.scalar)

  def test_product_sym_sym(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.SW1_full, self.SS_full)), self.TSW1 * self.TSS)

  def test_product_sym_full(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.SW1_full, self.R)), self.TSW1 * self.TR)

  def test_product_sym_skew(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.SW1_full, self.SW2_full)), self.TSW1 * self.TSW2)

  def test_product_skew_sym(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.SW1_full, self.WS_full)), self.TSW1 * self.TWS)

  def test_product_symmetric(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.SW1_full, self.S)), self.TSW1 * self.TS)

  def test_product_general(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.SW1_full, self.G)), self.TSW1 * self.TG)

  def test_product_skew(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.SW1_full, self.W)), self.TSW1 * self.TW)

class TestSkewSymR4(unittest.TestCase):
  def setUp(self):
    self.WS1 = np.array([
      [-8.3567359 , -5.39728818, -8.00844442, -8.33365112, -0.97903364, -8.23943149],
      [-6.97125417,  4.34802055,  7.06281056, -1.57511617,  7.83359933, -9.37625432],
      [-6.0799489 , -6.0309543 ,  3.68575895,  8.84296976,  6.55799427, -9.22029379]])
    self.WS1_full = common.wws2ts(self.WS1)
    self.TWS1 = tensors.SkewSymR4(self.WS1)

    self.WS2 = np.array([
      [-8.80662663,  0.46179936, -5.49454144,  7.91618428,  5.34053953, -6.68997405],
      [ 4.15874971, -4.59781751,  7.43746813,  8.99981425, -0.97692573, 2.5075246 ],
      [ 9.53201007, -8.03524224,  0.94329443, -6.44415877, -9.92911741, 3.51742689]])
    self.WS2_full = common.wws2ts(self.WS2)
    self.TWS2 = tensors.SkewSymR4(self.WS2)

    self.SW = np.array([
      [ 5.43434005, -6.55983214,  0.29737664],
      [-4.77472172, -8.51287287, -3.19380185],
      [ 4.43407952, -6.02555614,  5.87786914],
      [ 1.89488869, -5.65383917,  8.83717547],
      [-7.18030867,  1.56100537, -9.83238641],
      [-4.52369317, -3.07284914, -7.54966999]])
    self.SW_full = common.ws2ts(self.SW)
    self.TSW = tensors.SymSkewR4(self.SW)

    self.SS = np.array([
      [ 5.99159801, -2.24342348,  0.26667281, -0.95466199,  3.98931478, -0.10846981],
      [ 1.86468226, -4.32391908, -7.82738638, -7.45008989,  5.89874777, 0.45820648],
      [-5.92565398,  2.4862829 , -6.02112389,  6.75455965,  4.65183463, 9.96900579],
      [ 0.60378883, -3.72189328, -7.63388446, -5.76559403, -0.3119789 , -1.1527258 ],
      [ 4.56813135, -6.06783828, -6.18341368,  8.06169686, -9.56928844, 9.08114655],
      [-8.25516614,  6.30663846,  7.2084381 , -7.38280703, -5.96279902, 8.9935982 ]])
    self.SS_full = common.ms2ts(self.SS)
    self.TSS = tensors.SymSymR4(self.SS)

    self.R = np.array([[[[-8.03675620e+00,  2.58575052e+00,  2.44069661e+00],
             [ 4.75021663e+00,  1.24463394e+00, -8.69751301e-01],
             [-1.46310894e+00, -1.15053235e+00, -3.75342982e+00]],
            [[-7.64033956e+00,  4.19956720e+00, -4.87644982e+00],
             [ 1.06577507e+00,  8.94272637e+00,  6.57264250e-01],
             [-4.22613258e+00, -5.08830314e+00,  1.57718186e+00]],
            [[-4.02243082e+00, -4.75463781e+00, -8.88662152e+00],
             [-1.30383950e+00, -1.98063574e+00, -3.18963544e+00],
             [-7.52071674e+00,  1.08931933e+00,  2.86988431e+00]]],
           [[[ 5.28621060e+00, -6.83799668e+00,  8.98005935e+00],
             [-7.92741122e+00,  5.75699425e-01,  1.66782544e+00],
             [ 2.60041984e+00, -1.04476986e-02, -6.12424787e+00]],
            [[-3.73727368e+00,  6.59764771e+00, -1.18045587e+00],
             [ 4.08567441e+00,  2.66148943e+00, -6.82495588e-01],
             [-1.64417262e+00,  5.33119298e+00,  8.11045988e-03]],
            [[-5.90193883e+00, -2.63316107e+00,  5.61381825e+00],
             [-6.08591194e+00,  8.77285539e+00, -7.15230533e+00],
             [ 3.15093096e+00,  1.41350149e+00,  1.11702016e+00]]],
           [[[-9.61472764e-01, -1.91492497e+00,  9.48275324e+00],
             [ 6.68841134e+00,  3.23412041e+00, -3.41944541e+00],
             [-9.80203467e+00,  6.58425335e+00, -2.16548636e+00]],
            [[ 6.63950740e+00,  3.91551441e+00, -8.98229111e+00],
             [ 9.84606756e+00, -8.16145090e+00,  8.41929062e-01],
             [-1.93839620e+00,  7.44485127e+00, -2.70832414e+00]],
            [[ 9.79265531e+00, -1.18212395e+00, -5.39433704e+00],
             [ 4.87152614e+00,  9.47287450e+00,  5.53838514e+00],
             [ 9.30443367e+00,  1.27090319e+00,  1.60409739e+00]]]])
    self.TR = tensors.RankFour(self.R)

    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)

    self.scalar = 5.2

    self.G = np.array([[ 9.50640677,  1.79084726, -2.8877036 ],
       [-1.63159958,  2.52866904, -8.71585042],
       [ 5.01859685, -8.7324075 , -0.42919134]])
    self.TG = tensors.RankTwo(self.G)

    self.W = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.W = 0.5 * (self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

  def test_to_full(self):
    full_np = common.wws2ts(self.WS1)
    full_t = tensors.RankFour(full_np)
    full = self.TWS1.to_full()
    self.assertEqual(full_t, full)

  def test_from_full(self):
    full = self.TWS1.to_full()
    new = full.to_skewsym()
    self.assertEqual(self.TWS1, new)

  def test_add(self):
    self.assertEqual(tensors.SkewSymR4(self.WS1 + self.WS2), self.TWS2 + self.TWS1)
    self.assertEqual(tensors.SkewSymR4(self.WS1 - self.WS2), self.TWS1 - self.TWS2)

  def test_equality(self):
    self.assertEqual(self.TWS1, self.TWS1)

  def test_inequality(self):
    self.assertNotEqual(self.TWS1, self.TWS2)

  def test_negate(self):
    self.assertEqual(tensors.SkewSymR4(-self.WS1), -self.TWS1)

  def test_scalar_mult(self):
    self.assertEqual(tensors.SkewSymR4(self.scalar * self.WS1), self.scalar * self.TWS1)
    self.assertEqual(tensors.SkewSymR4(self.scalar * self.WS2), self.TWS2 * self.scalar)
    self.assertEqual(tensors.SkewSymR4(self.WS1 / self.scalar), self.TWS1 / self.scalar)

  def test_product_sym_sym(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.WS1_full, self.SS_full)), self.TWS1 * self.TSS)

  def test_product_sym_full(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.WS1_full, self.R)), self.TWS1 * self.TR)

  def test_product_sym_skew(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.WS1_full, self.SW_full)), self.TWS1 * self.TSW)

  def test_product_skew_sym(self):
    self.assertEqual(tensors.RankFour(np.einsum('ijkl,klmn', self.WS1_full, self.WS2_full)), self.TWS1 * self.TWS2)

  def test_product_symmetric(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.WS1_full, self.S)), self.TWS1 * self.TS)

  def test_product_general(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.WS1_full, self.G)), self.TWS1 * self.TG)

  def test_product_skew(self):
    self.assertEqual(tensors.RankTwo(np.einsum('ijkl,kl', self.WS1_full, self.W)), self.TWS1 * self.TW)

  def test_douter(self):
    self.assertEqual(tensors.SkewSymR4(common.ts2wws(np.einsum('ij,kl', self.W, self.S))), tensors.douter(self.TW, self.TS))

class TestCPSpeciality(unittest.TestCase):
  def setUp(self):
    self.SS = np.array([
      [ 5.99159801, -2.24342348,  0.26667281, -0.95466199,  3.98931478, -0.10846981],
      [ 1.86468226, -4.32391908, -7.82738638, -7.45008989,  5.89874777, 0.45820648],
      [-5.92565398,  2.4862829 , -6.02112389,  6.75455965,  4.65183463, 9.96900579],
      [ 0.60378883, -3.72189328, -7.63388446, -5.76559403, -0.3119789 , -1.1527258 ],
      [ 4.56813135, -6.06783828, -6.18341368,  8.06169686, -9.56928844, 9.08114655],
      [-8.25516614,  6.30663846,  7.2084381 , -7.38280703, -5.96279902, 8.9935982 ]])
    self.SS_full = common.ms2ts(self.SS)
    self.TSS = tensors.SymSymR4(self.SS)

    self.W = np.array([[-9.36416517,  2.95527444,  8.70983194],
           [-1.54693052,  8.7905658 , -5.10895168],
           [-8.52740468, -0.7741642 ,  2.89544992]])
    self.W = 0.5 * (self.W - self.W.T)
    self.TW = tensors.Skew(self.W)

    self.S = np.array([[4.1,2.8,-1.2],[3.1,7.1,0.2],[4,2,3]])
    self.S = 0.5*(self.S + self.S.T)
    self.TS = tensors.Symmetric(self.S)

    self.WS = np.array([
      [-8.3567359 , -5.39728818, -8.00844442, -8.33365112, -0.97903364, -8.23943149],
      [-6.97125417,  4.34802055,  7.06281056, -1.57511617,  7.83359933, -9.37625432],
      [-6.0799489 , -6.0309543 ,  3.68575895,  8.84296976,  6.55799427, -9.22029379]])
    self.WS_full = common.wws2ts(self.WS)
    self.TWS = tensors.SkewSymR4(self.WS)

  def test_symsymskew_skewsymsym(self):
    A1 = tensors.SymSymR4Skew_SkewSymR4SymR4(self.TSS, self.TW)

    A2_ten = np.einsum('kmst,ml', self.SS_full, self.W) - np.einsum('km,mlst',
        self.W, self.SS_full)
    A2 = tensors.SymSymR4(common.ts2ms(A2_ten))

    self.assertEqual(A1, A2)

  def test_symskewsym_skewsymsym(self):
    A1 = tensors.SymSkewR4Sym_SkewSymR4SymR4(self.TWS, self.TS)

    A2_ten = np.einsum('km,mlst', self.S, self.WS_full) - np.einsum('kmst,ml',
        self.WS_full, self.S)
    A2 = tensors.SymSymR4(common.ts2ms(A2_ten))

    self.assertEqual(A1, A2)

  def test_special(self):
    A1 = tensors.SpecialSymSymR4Sym(self.TSS, self.TS)
    A2_ten = np.einsum('ijkz,ky', self.SS_full, self.S) - np.einsum(
        'ijyl,zl', self.SS_full, self.S)
    A2 = tensors.SymSkewR4(common.ts2sww(A2_ten))

    self.assertEqual(A1, A2)

class TestSymSymSymR6(unittest.TestCase):
  def setUp(self):
    self.SS = np.array([
      [ 5.99159801, -2.24342348,  0.26667281, -0.95466199,  3.98931478, -0.10846981],
      [ 1.86468226, -4.32391908, -7.82738638, -7.45008989,  5.89874777, 0.45820648],
      [-5.92565398,  2.4862829 , -6.02112389,  6.75455965,  4.65183463, 9.96900579],
      [ 0.60378883, -3.72189328, -7.63388446, -5.76559403, -0.3119789 , -1.1527258 ],
      [ 4.56813135, -6.06783828, -6.18341368,  8.06169686, -9.56928844, 9.08114655],
      [-8.25516614,  6.30663846,  7.2084381 , -7.38280703, -5.96279902, 8.9935982 ]])
    self.SS_full = common.ms2ts(self.SS)
    self.TSS = tensors.SymSymR4(self.SS)

    self.S = np.array([[ 6.19999242, -6.95811611, -6.02901899],
           [ 8.38508084,  6.01607694,  6.79839425],
           [-4.4214246 , -2.36795313, -8.84070728]])
    self.S = 0.5*(self.S+self.S.T)
    self.TS = tensors.Symmetric(self.S)

    self.FS = ra.random((6,6,6))
    self.FS_full = common.m62t6(self.FS)
    self.TFS = tensors.SymSymSymR6(self.FS)

  def test_dot_i(self):
    full = np.einsum('ijklmn,ij', self.FS_full, self.S)
    vec = self.TFS.dot_i(self.TS)
    vec2f = common.ms2ts(vec.data.reshape(6,6))

    self.assertTrue(np.allclose(full,vec2f))

  def test_dot_j(self):
    full = np.einsum('ijklmn,kl', self.FS_full, self.S)
    vec = self.TFS.dot_j(self.TS)
    vec2f = common.ms2ts(vec.data.reshape(6,6))

    self.assertTrue(np.allclose(full,vec2f))

  def test_dot_k(self):
    full = np.einsum('ijklmn,mn', self.FS_full, self.S)
    vec = self.TFS.dot_k(self.TS)
    vec2f = common.ms2ts(vec.data.reshape(6,6))

    self.assertTrue(np.allclose(full,vec2f))

  def test_outer_project(self):
    full = np.einsum('ijkl,mn', self.SS_full, self.S)

    vec = tensors.outer_product_k(self.TSS, self.TS)
    vec2f = common.m62t6(vec.data.reshape((6,6,6)))
    
    self.assertTrue(np.allclose(full, vec2f))

  def test_add(self):
    self.assertTrue(np.allclose((self.TFS+self.TFS).data, (2*self.FS).flatten()))

  def test_subtract(self):
    self.assertTrue(np.allclose((self.TFS-self.TFS).data, np.zeros((6,6,6)).flatten()))

  def test_negate(self):
    self.assertTrue(np.allclose((-self.TFS).data, -self.FS.flatten()))

  def test_scalar_mult(self):
    self.assertTrue(np.allclose((2.0*self.TFS).data, 2*self.FS.flatten()))
    self.assertTrue(np.allclose((self.TFS*2.0).data, 2*self.FS.flatten()))

  def test_scalar_divide(self):
    self.assertTrue(np.allclose((self.TFS/2.0).data, self.FS.flatten()/2))

