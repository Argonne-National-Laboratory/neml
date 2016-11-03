import sys
sys.path.append('..')

from neml import neml, volumetric, deviatoric, shear

from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonMatModel(object):
  """
    Tests that could be applied to all material models
  """
  def test_history(self):
    self.assertEqual(self.model.nstore, len(self.hist0))

  def test_tangent_sd(self):
    s_np1, h_np1, A_np1 = self.model.update_sd(self.e_np1, self.e_n,
        self.T_np1, self.T_n,
        self.t_np1, self.t_n,
        self.s_n, self.h_n)
    dfn = lambda x: self.model.update_sd(x, self.e_n,
        self.T_np1, self.T_n,
        self.t_np1, self.t_n,
        self.s_n, self.h_n)[0]
    n_A = differentiate(dfn, self.e_np1)
    self.assertTrue(np.allclose(n_A, A_np1),
        msg = str(A_np1))

class TestLinearElasticSplit(unittest.TestCase, CommonMatModel):
  def setUp(self):
    self.E = 120000.0
    self.nu = 0.3

    self.mu = self.E / (2 * (1 + self.nu))
    self.K = self.E / (3 * (1 - 2 * self.nu))

    vol_model = volumetric.VModelConstantK(self.K)
    shear_model = shear.ConstantShearModulus(self.mu)
    dev_model = deviatoric.LEModel(shear_model)
    self.model = neml.SplitModel_sd(vol_model, dev_model)

    self.hist0 = np.array([])
    self.h_n = np.array([])

    self.s_n = np.array([50.0,-20.0,30.0,10.0,5.0,-15.0])
    self.e_n = np.array([0.1,-0.05,-0.025,0.01,0.02,-0.0025])
    self.e_np1 = np.array([0.09,0.04,0.04,0.02,0.02,0.1])

    self.T_np1 = 300.0
    self.T_n = 270.0

    self.t_np1 = 1.0
    self.t_n = 0.0

    self.C = self.E / ((1+self.nu) * (1-2*self.nu)) * np.array([
      [1-self.nu, self.nu, self.nu, 0, 0 ,0],
      [self.nu, 1-self.nu, self.nu, 0, 0, 0],
      [self.nu, self.nu, 1-self.nu, 0, 0, 0],
      [0, 0, 0, (1-2.0*self.nu), 0, 0],
      [0, 0, 0, 0, (1-2.0*self.nu), 0],
      [0, 0, 0, 0, 0, (1-2.0*self.nu)]])

  def test_update(self):
    s_np1, h_np1, A_np1 = self.model.update_sd(self.e_np1, self.e_n,
        self.T_np1, self.T_n,
        self.t_np1, self.t_n,
        self.s_n, self.h_n)
    self.assertTrue(np.allclose(s_np1, np.dot(self.C, self.e_np1)),
        msg = str(np.dot(self.C, self.e_np1)))
