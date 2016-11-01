import sys
sys.path.append('..')

from neml.nemlmath import *
from common import *

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

class TestInvert(unittest.TestCase):
  def setUp(self):
    self.n = 10
    self.square_ns = ra.random((self.n,self.n))
    self.nonsquare = ra.random((self.n,self.n-1))
    self.big = ra.random((self.n,self.n,self.n))

  def test_invert(self):
    inv_ns = la.inv(self.square_ns)
    inv = invert_matrix(self.square_ns)
    self.assertTrue(np.allclose(inv_ns, inv))

  def test_nonsquare(self):
    self.assertRaises(RuntimeError, invert_matrix, self.nonsquare)

  def test_nonmatrix(self):
    self.assertRaises(RuntimeError, invert_matrix, self.big)
