import sys
sys.path.append('..')

from neml import deviatoric, shear, surface, hardening

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

class CommonPlasticity(object):
  """
    Common testing protocol for history-dependent models.
    You need to set e_max, t_max, and nsteps
    You want e_max to get well into the plastic range and nsteps enough
    to keep the model convergent
  """
  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_tangent_steps(self):
    """
      Tangent test
    """
    h_n = self.model.init_hist()
    s_n = np.zeros((6,))
    e_n = np.zeros((6,))
    t_n = 0.0
    for i,f in enumerate(np.linspace(0,1,self.nsteps)):
      print(i)
      e_np1 = self.e_max * f
      t_np1 = self.t_max * f

      s_np1, h_np1, A_np1 = self.model.update(e_np1, e_n,
          self.T, self.T, t_np1, t_n, s_n, h_n)

      dfn = lambda x: self.model.update(x, e_n, 
          self.T, self.T, t_np1, t_n, s_n, h_n)[0]
      num_A = differentiate(dfn, e_np1)
      
      self.assertTrue(np.allclose(A_np1, num_A, rtol = 1.0e-6),
          msg = str((A_np1)))

      e_n = np.copy(e_np1)
      t_n = t_np1
      s_n = np.copy(s_np1)
      h_n = np.copy(h_np1)

    #self.assertTrue(False)


class AssociativeLinearKin(unittest.TestCase, CommonPlasticity):
  def setUp(self):
    self.e_max = np.array([0.05,0,0,0,0,0])/10.0
    self.t_max = 10.0
    self.nsteps = 100

    self.T = 300.0

    self.E = 200000.0
    self.nu = 0.3

    self.mu = self.E / (2 * (1.0 + self.nu))
    self.K = self.E / (3 * (1 - 2 * self.nu))
    
    self.K0 = 200.0
    self.Kp = self.E / 100.0

    shear_model = shear.ConstantShearModulus(self.mu)
    self.ys = surface.IsoJ2()
    hr = hardening.IsoJ2LinearAHardening(self.K0, self.Kp)
    self.model = deviatoric.RIAFModel(shear_model, self.ys, hr, verbose=True)

    self.hist0 = np.zeros((7,))



