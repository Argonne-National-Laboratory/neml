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

class TestOuter(unittest.TestCase):
  def setUp(self):
    self.na = 10
    self.nb = 12
    self.a = ra.random((self.na,))
    self.b = ra.random((self.nb,))

  def test_outer(self):
    self.assertTrue(
        np.allclose(outer_vec(self.a, self.b), np.outer(self.a, self.b)),
        msg = str(outer_vec(self.a, self.b)))

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
