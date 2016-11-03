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
    self.e_np1 = np.array([0.01,-0.02,0.03,0.1,-0.2,0.01])
    self.e_n = np.array([0.0,0.001,-0.01,0.11,0.05,0.1])
    self.e_np1_mean = np.sum(self.e_np1[:3])
    self.e_n_mean = np.sum(self.e_n[:3])
    self.h_n = np.array([])
    self.s_n = np.zeros((6,))
    self.s_n_mean = np.sum(self.s_n[:3])
    
    self.T_np1 = 300.0
    self.T_n = 270

    self.t_np1 = 1.0
    self.t_n = 0.1

  def test_properties(self):
    self.assertTrue(np.isclose(self.K, self.model.K))

  def test_history(self):
    self.assertEqual(self.model.nhist, 0)
    self.assertTrue(np.allclose(self.model.init_hist(), np.array([])))

  def test_mean_update(self):
    s_np1, h_np1, A_np1 = self.model.update_mean(self.e_np1_mean, self.e_n_mean,
        self.T_np1, self.T_n, self.t_np1, self.t_n,
        self.s_n_mean, self.h_n)
    self.assertTrue(np.isclose(s_np1, self.K * self.e_np1_mean))

  def test_mean_update_derivative(self):
    s_np1, h_np1, A_np1 = self.model.update_mean(self.e_np1_mean, self.e_n_mean,
        self.T_np1, self.T_n, self.t_np1, self.t_n,
        self.s_n_mean, self.h_n)

    dfn = lambda x: self.model.update_mean(x, self.e_n_mean, self.T_np1, self.T_n,
        self.t_np1, self.t_n, self.s_n_mean, self.h_n)[0]

    n_A = differentiate(dfn, self.e_np1_mean)

    self.assertTrue(np.isclose(n_A, A_np1))

  def test_update(self):
    s_np1, h_np1, A_np1 = self.model.update(self.e_np1, self.e_n,
        self.T_np1, self.T_n, self.t_np1, self.t_n,
        self.s_n, self.h_n)

    self.assertTrue(np.allclose(s_np1, np.array([1,1,1,0,0,0]) * self.K * 
      self.e_np1_mean))

  def test_update_derivative(self):
    s_np1, h_np1, A_np1 = self.model.update(self.e_np1, self.e_n,
        self.T_np1, self.T_n, self.t_np1, self.t_n,
        self.s_n, self.h_n)

    dfn = lambda x: self.model.update(x, self.e_n, self.T_np1, self.T_n,
        self.t_np1, self.t_n, self.s_n, self.h_n)[0]

    n_a = differentiate(dfn, self.e_np1)

    self.assertTrue(np.allclose(n_a, A_np1))
