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


