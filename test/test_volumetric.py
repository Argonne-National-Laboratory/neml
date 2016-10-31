import sys
sys.path.append('..')

from neml.volumetric import *
from common import *

import unittest
import numpy as np

class TestVModelConstantK(unittest.TestCase):
  def setUp(self):
    self.K = 1000.0
    self.model = VModelConstantK(self.K)
    self.e_inc = np.array([0.01,-0.02,0.03,0.1,-0.2,0.01])
    self.e_inc_mean = np.sum(self.e_inc[:3])
    self.h_n = np.array([])
    self.s_n = np.zeros((6,))
    self.s_n_mean = np.sum(self.s_n[:3])
    
    self.T_np1 = 300.0
    self.T_inc = 30.0

    self.t_np1 = 0.0
    self.t_inc = 0.1

  def test_properties(self):
    self.assertTrue(np.isclose(self.K, self.model.K))

  def test_mean_update(self):
    h_np1, s_np1, A_np1 = self.model.update_mean(self.e_inc_mean,
        self.T_np1, self.T_inc, self.t_np1, self.t_inc,
        self.h_n, self.s_n_mean)
    self.assertTrue(np.isclose(s_np1, self.s_n_mean + self.K * self.e_inc_mean))

  def test_mean_update_derivative(self):
    h_np1, s_np1, A_np1 = self.model.update_mean(self.e_inc_mean,
        self.T_np1, self.T_inc, self.t_np1, self.t_inc,
        self.h_n, self.s_n_mean)

    dfn = lambda x: self.model.update_mean(x, self.T_np1, self.T_inc,
        self.t_np1, self.t_inc, self.h_n, self.s_n_mean)[1]

    n_A = differentiate(dfn, self.e_inc_mean)

    self.assertTrue(np.isclose(n_A, A_np1))

  def test_update(self):
    h_np1, s_np1, A_np1 = self.model.update(self.e_inc,
        self.T_np1, self.T_inc, self.t_np1, self.t_inc,
        self.h_n, self.s_n)

    self.assertTrue(np.allclose(s_np1, np.array([1,1,1,0,0,0]) * self.K * 
      self.e_inc_mean))

  def test_update_derivative(self):
    h_np1, s_np1, A_np1 = self.model.update(self.e_inc,
        self.T_np1, self.T_inc, self.t_np1, self.t_inc,
        self.h_n, self.s_n)

    dfn = lambda x: self.model.update(x, self.T_np1, self.T_inc,
        self.t_np1, self.t_inc, self.h_n, self.s_n)[1]

    n_a = differentiate(dfn, self.e_inc)

    self.assertTrue(np.allclose(n_a, A_np1))
