import sys
sys.path.append('..')

from neml import deviatoric, shear

from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonDeviatoric(object):
  """
    Common deviatoric tests
  """
  def test_history(self):
    self.assertEqual(len(self.hist0), self.model.nhist)
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_tangent(self):
    s_np1, h_np1, A_np1 = self.model.update(self.e_np1, self.e_n,
        self.T_np1, self.T_n,
        self.t_np1, self.t_n,
        self.s_n, self.h_n)
    dfn = lambda x: self.model.update(x, self.e_n,
        self.T_np1, self.T_n,
        self.t_np1, self.t_n,
        self.s_n, self.h_n)[0]
    n_A = differentiate(dfn, self.e_np1)
    self.assertTrue(np.allclose(n_A, A_np1),
        msg = str(A_np1))

class LEModel(unittest.TestCase, CommonDeviatoric):
  def setUp(self):
    self.s_n = np.array([50.0,-20.0,30.0,10.0,5.0,-15.0])
    self.e_n = np.array([0.1,-0.05,-0.025,0.01,0.02,-0.0025])
    self.e_np1 = np.array([0.09,0.04,0.0,0.02,0.0,0.1])

    self.e_np1_dev = self.e_np1 - np.array([1,1,1,0,0,0]) * np.sum(self.e_np1[:3]) / 3.0

    self.s_n = np.array([50,20.0,10.0,5.0,2.0,-10.0])

    self.T_np1 = 300.0
    self.T_n = 270.0

    self.t_np1 = 1.0
    self.t_n = 0.0

    self.hist0 = np.array([])
    self.h_n = np.array([])

    self.mu = 120000.0
    smodel = shear.ConstantShearModulus(self.mu)
    self.model = deviatoric.LEModel(smodel)

  def test_update(self):
    s_np1, h_np1, A_np1 = self.model.update(self.e_np1, self.e_n,
        self.T_np1, self.T_n,
        self.t_np1, self.t_n,
        self.s_n, self.h_n)
    self.assertTrue(np.allclose(s_np1, 2.0 * self.mu * self.e_np1_dev))

