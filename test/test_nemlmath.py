import sys
sys.path.append('..')

from neml.nemlmath import *
from common import *

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

class TestAddSubtractNegate(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.a = ra.random((self.n,))
    self.b = ra.random((self.n,))

  def test_minus(self):
    acopy = np.copy(self.a)
    self.assertTrue(np.allclose(minus_vec(self.a), -acopy))

  def test_add(self):
    self.assertTrue(np.allclose(add_vec(self.a, self.b), self.a + self.b))

  def test_sub(self):
    self.assertTrue(np.allclose(sub_vec(self.a, self.b), self.a - self.b))

class TestDotNorm(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.a = ra.random((self.n,))
    self.b = ra.random((self.n,))

  def test_dot(self):
    self.assertTrue(np.isclose(dot_vec(self.a,self.b), np.dot(self.a,self.b)))

  def test_norm(self):
    self.assertTrue(np.isclose(norm2_vec(self.a), la.norm(self.a)))

  def test_normalize(self):
    nv = self.a / la.norm(self.a)
    self.assertTrue(np.allclose(nv, normalize_vec(self.a)))

class TestDev(unittest.TestCase):
  def setUp(self):
    self.s = np.array([10.0,-15.0,5.0,4.0,8.0,-12.0])

  def test_dev(self):
    sdev = self.s - np.array([1,1,1,0,0,0]) * np.sum(self.s[:3]) / 3.0
    self.assertTrue(np.allclose(sdev, dev_vec(self.s)))

class TestOuter(unittest.TestCase):
  def setUp(self):
    self.na = 10
    self.nb = 12
    self.a = ra.random((self.na,))
    self.b = ra.random((self.nb,))
    self.C = ra.random((self.na,self.nb))

  def test_outer(self):
    self.assertTrue(
        np.allclose(outer_vec(self.a, self.b), np.outer(self.a, self.b)),
        msg = str(outer_vec(self.a, self.b)))

  def test_update(self):
    Cu = self.C + np.outer(self.a, self.b)
    self.assertTrue(
        np.allclose(outer_update(self.a, self.b, self.C), Cu))

  def test_update_minus(self):
    Cu = self.C - np.outer(self.a, self.b)
    self.assertTrue(
        np.allclose(outer_update_minus(self.a, self.b, self.C), Cu))

class TestMatVec(unittest.TestCase):
  def setUp(self):
    self.m = 10
    self.n = 8
    self.A = ra.random((self.m, self.n))
    self.b = ra.random((self.n,))
    self.c = ra.random((self.m,))

  def test_matvec(self):
    self.assertTrue(np.allclose(np.dot(self.A, self.b), 
      mat_vec(self.A, self.b)),
      msg = str(np.dot(self.A, self.b)) + '\n' + str(mat_vec(self.A, self.b)))

  def test_matvec_trans(self):
    self.assertTrue(np.allclose(np.dot(self.A.T, self.c),
      mat_vec_trans(self.A, self.c)))

class TestMatMat(unittest.TestCase):
  def setUp(self):
    self.m = 10
    self.n = 8
    self.k = 9
    self.A = ra.random((self.m, self.k))
    self.B = ra.random((self.k, self.n))

  def test_matmat(self):
    self.assertTrue(np.allclose(np.dot(self.A, self.B), mat_mat(self.A, self.B)))

class TestInvert(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.square_ns = ra.random((self.n,self.n))
    self.nonsquare = ra.random((self.n,self.n-1))
    self.big = ra.random((self.n,self.n,self.n))

  def test_invert(self):
    inv_ns = la.inv(self.square_ns)
    inv = invert_mat(self.square_ns)
    self.assertTrue(np.allclose(inv_ns, inv))

  def test_nonsquare(self):
    self.assertRaises(RuntimeError, invert_mat, self.nonsquare)

  def test_nonmatrix(self):
    self.assertRaises(RuntimeError, invert_mat, self.big)

class TestSolve(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.A = ra.random((self.n,self.n))
    self.b = ra.random((self.n,))

  def test_solve(self):
    x = la.solve(self.A, self.b)
    self.b = solve_mat(self.A, self.b)
    print(x)
    print(self.b)
    self.assertTrue(np.allclose(x, self.b))

class TestCond(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.A = ra.random((self.n,self.n))

  def test_cond(self):
    fr = condition(self.A)
