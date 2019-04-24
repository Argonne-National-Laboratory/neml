import sys
sys.path.append('..')

from neml.nemlmath import *
from common import *

from neml import nemlmath

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

class TestLDTangent(unittest.TestCase):
  def setUp(self):
    self.sym_mandel = ra.random((6,6))
    self.sym_ten = ms2ts(self.sym_mandel)
    self.sym_99 = unroll_fourth(self.sym_ten)

    self.skew_deriv_red = ra.random((6,3))
    self.skew_deriv_ten = ws2ts(self.skew_deriv_red)
    self.skew_deriv_99 = unroll_fourth(self.skew_deriv_ten)

    self.S = np.array([[10.0,-20.0,15.0],[-20.0,60.0,5.0],[15.0,5.0,20.0]])
    self.S_v = sym(self.S)

  def test_piece_back(self):
    A = transform_fourth(self.sym_mandel, self.skew_deriv_red).reshape((3,3,3,3))
    B = piece_together_fourth(self.sym_mandel, self.skew_deriv_red)

    self.assertTrue(np.allclose(A, B))

  def test_mandel_convert(self):
    self.assertTrue(np.allclose(nemlmath.mandel2full(self.sym_mandel), self.sym_99))
    self.assertTrue(np.allclose(nemlmath.full2mandel(self.sym_99), self.sym_mandel))

  def test_skew_convert(self):
    self.assertTrue(np.allclose(nemlmath.skew2full(self.skew_deriv_red), self.skew_deriv_99))
    self.assertTrue(np.allclose(nemlmath.full2skew(self.skew_deriv_99), self.skew_deriv_red))

  def test_identities(self):
    sym = unroll_fourth(0.5*(np.einsum('ik,jl', np.eye(3), np.eye(3)) + np.einsum('jk,il', np.eye(3), np.eye(3))))
    skew = unroll_fourth(0.5*(np.einsum('ik,jl', np.eye(3), np.eye(3)) - np.einsum('jk,il', np.eye(3), np.eye(3)))) 

    self.assertTrue(np.allclose(sym, nemlmath.idsym()))
    self.assertTrue(np.allclose(skew, nemlmath.idskew()))

  def test_outer_mat(self):
    tensor = np.einsum('in,jm', self.S, np.eye(3)) + np.einsum('im,nj', np.eye(3), self.S) - np.einsum('mn,ij', np.eye(3), self.S)
    matrix = unroll_fourth(tensor)
    
    direct = nemlmath.truesdell_tangent_outer(self.S_v)
    self.assertTrue(np.allclose(matrix, direct))

class TestLDUpdate(unittest.TestCase):
  def setUp(self):
    self.L = np.array([[0.1,-0.2,0.5],[1.1,-0.2,0.25],[0.3,-0.4,0.2]])
    self.D = 0.5 * (self.L + self.L.T)
    self.Dv = sym(self.D)
    self.W = 0.5 * (self.L - self.L.T)
    self.Wv = skew(self.W)
    
    self.S_n = np.array([[10.0,-20.0,15.0],[-20.0,60.0,5.0],[15.0,5.0,20.0]])
    self.S_n_v = sym(self.S_n)

    self.dS = np.array([[30.0,10.0,5.0],[10.0,-15.0,15.0],[5.0,15.0,-10.0]])
    self.dS_v = sym(self.dS)

  def test_setup_right(self):
    self.assertTrue(np.allclose(self.S_n, self.S_n.T))
    self.assertTrue(np.allclose(self.dS, self.dS.T))

  def test_gives_sym(self):
    """
      Make sure the update actually gives a symmetric tensor
    """
    val = self.dS - np.einsum('ik,jk', self.S_n, self.L) - np.einsum('ik,kj', self.L, self.S_n) + np.trace(self.L
        ) * self.S_n
    self.assertTrue(np.allclose(val, val.T)) 

  def test_update_theory(self):
    """
      Test to make sure my tensor math was correct on deriving the form of the update (tensor notation)
    """
    J = np.einsum('im,jn', np.eye(3), np.eye(3))*(1.0 + np.trace(self.L)) - np.einsum('im,jn', np.eye(3), self.L) - np.einsum('im,jn',
        self.L, np.eye(3))
    Jinv = la.tensorinv(J)
    interm = self.dS + np.dot(self.S_n, self.L.T) + np.dot(self.L, self.S_n) - np.trace(self.L) * self.S_n
    S_np1 = self.S_n + np.einsum('ijkl,kl', Jinv, interm)
    
    dS_c = S_np1 - self.S_n

    verify = dS_c - np.dot(S_np1, self.L.T) - np.dot(self.L, S_np1) + np.trace(self.L) * S_np1
    self.assertTrue(np.allclose(self.dS, verify))

  def test_to_vector(self):
    """
      Make sure the C++ tensor -> Mandel works
    """
    self.assertTrue(np.allclose(self.Dv, nemlmath.sym(self.D)))

  def test_to_tensor(self):
    """
      Make sure the C++ Mandel -> tensor works
    """
    self.assertTrue(np.allclose(self.D, nemlmath.usym(self.Dv)))

  def test_rhs(self):
    """
      Make sure the Truesdell rhs is correct
    """
    interm = sym(self.dS + np.dot(self.S_n, self.L.T) + np.dot(self.L, self.S_n) - np.trace(self.L) * self.S_n)
    outterm = truesdell_rhs(self.Dv, self.Wv, self.S_n_v, self.dS_v)
    self.assertTrue(np.allclose(interm, outterm))

  def test_matrix(self):
    """
      Make sure the Truesdell matrix is correct
    """
    J = np.einsum('im,jn', np.eye(3), np.eye(3))*(1.0 + np.trace(self.L)) - np.einsum('im,jn', np.eye(3), self.L) - np.einsum('im,jn',
        self.L, np.eye(3))
    Jm = unroll_fourth(J)
    
    j2 = truesdell_mat(self.Dv, self.Wv)

    self.assertTrue(np.allclose(Jm, j2))

  def test_implemented_update(self):
    """
      Test to make sure the vector update function will work
    """
    J = np.einsum('im,jn', np.eye(3), np.eye(3))*(1.0 + np.trace(self.L)) - np.einsum('im,jn', np.eye(3), self.L) - np.einsum('im,jn',
        self.L, np.eye(3))
    Jinv = la.tensorinv(J)
    interm = self.dS + np.dot(self.S_n, self.L.T) + np.dot(self.L, self.S_n) - np.trace(self.L) * self.S_n
    S_np1 = self.S_n + np.einsum('ijkl,kl', Jinv, interm)
    
    S_np1_v = sym(S_np1)

    S_np1_2 = truesdell_update_sym(self.Dv, self.Wv, self.S_n_v, self.dS_v)

    self.assertTrue(np.allclose(S_np1_v, S_np1_2))

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
