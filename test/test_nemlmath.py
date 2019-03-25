import sys
sys.path.append('..')

from neml.nemlmath import *
from common import *

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

class TestThingsWithSkew(unittest.TestCase):
  def setUp(self):
    self.B1 = np.array([[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0]])
    self.E = 0.5*(self.B1 + self.B1.T)
    self.Ev = sym(self.E)
    self.B2 = np.array([[-3.0,-1.0,2.0],[-4.0,5.0,6.0],[9.0,-7.0,-8.0]])
    self.W = 0.5*(self.B2 - self.B2.T)
    self.Wv = skew(self.W)

    self.E2 = 0.5*(self.B2 + self.B2.T)
    self.E2v = sym(self.E2)

    self.B3 = np.array([[0.1,0.5,0.2],[-0.25,0.6,0.4],[1.0,3.0,1.0]])
    self.E3 = 0.5*(self.B3 + self.B3.T)

  def test_my_understanding(self):
    ten = np.einsum('im,jn', np.eye(3), np.eye(3)) + np.einsum('im,nj', np.eye(3), self.W) - np.einsum('im,jn', self.W, np.eye(3))
    iten = la.tensorinv(ten)
    
    issym = np.einsum('ijkl,kl', iten, self.E2)
    self.assertTrue(np.allclose(issym - issym.T, np.zeros((3,3))))

    # None of our existing matrix routines will work
    # We will need a custom operator for the action of the inverse tensor on a Mandel vector

    # How about for general matrices?
    ten = np.einsum('im,jn', np.eye(3), np.eye(3)) + np.einsum('im,nj', np.eye(3), self.B1) - np.einsum('im,jn', self.B1, np.eye(3))
    iten = la.tensorinv(ten)
    
    issym = np.einsum('ijkl,kl', iten, self.E2)
    self.assertFalse(np.allclose(issym - issym.T, np.zeros((3,3))))

    # So what about the full Lie derivative?  Maybe that extra tr term is needed?

  def test_skew_update_formula(self):
    # Check how I plan do to the update (long way)
    W = self.W
    dS_jau = self.E
    Sn = self.E2

    J = np.einsum('im,jn', np.eye(3), np.eye(3)) + np.einsum('im,nj', np.eye(3), W) - np.einsum('im,jn', W, np.eye(3))
    Jinv = la.tensorinv(J)

    T1 = dS_jau - np.einsum('mk,kn', Sn, W) + np.einsum('mk,kn', W, Sn)

    Snp1 = Sn + np.einsum('ijmn,mn', Jinv, T1)

    dS1 = Snp1 - Sn

    dS2 = dS_jau - np.einsum('ik,kj', Snp1, W) + np.einsum('ik,kj', W, Snp1)

    self.assertTrue(np.allclose(dS1, dS2))

  def test_inverse_formula(self):
    """
      Check the direct inverse formula
    """
    W = self.W

    J = np.einsum('im,jn', np.eye(3), np.eye(3)) + np.einsum('im,nj', np.eye(3), W) - np.einsum('im,jn', W, np.eye(3))
    Jinv = la.tensorinv(J)

    print(Jinv)
    self.assertTrue(False)

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
    Cn = outer_update(self.a, self.b, self.C)
    self.assertTrue(np.allclose(Cn, Cu))

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

class TestDiagSolve(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.dl = ra.random((self.n-1,))
    self.d = ra.random((self.n,))
    self.du = ra.random((self.n-1))
    self.A = np.zeros((self.n, self.n))
    for i in range(self.n):
      self.A[i,i] = self.d[i]
      if i != self.n -1:
        self.A[i,i+1] = self.du[i]
        self.A[i+1,i] = self.dl[i]
    
    self.c = ra.random((self.n,))

  def test_solve(self):
    x_np = la.solve(self.A, self.c)
    dl, d, du, du2, ipiv = dgttrf(self.dl, self.d, self.du)
    x_la = dgttrs(dl, d, du, du2, ipiv, self.c)
    self.assertTrue(np.allclose(x_np, x_la))

class TestCond(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.A = ra.random((self.n,self.n))

  def test_cond(self):
    fr = condition(self.A)

class TestPoly(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.poly = ra.random((self.n,))
    self.x = ra.random((1,))[0]
  
  def test_eval(self):
    npv = np.polyval(self.poly, self.x)
    mv  = polyval(self.poly, self.x)
    self.assertTrue(np.isclose(npv, mv))
