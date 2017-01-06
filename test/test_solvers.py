import sys
sys.path.append('..')

from neml import solvers
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonSolvable(object):
  """
    Common tests for a Solvable object
  """
  def test_jacobian(self):
    x = self.gen_x()
    dummy = solvers.TrialState()
    R, J = self.model.RJ(x, dummy)
    dfn = lambda x: self.model.RJ(x, dummy)[0]
    nJ = differentiate(dfn, x)

    self.assertTrue(np.allclose(J, nJ, rtol = 1.0e-3, atol = 1.0e-4))

class TestRosenbrock(unittest.TestCase, CommonSolvable):
  def setUp(self):
    self.N = 5
    self.model = solvers.TestRosenbrock(self.N)

  def gen_x(self):
    return np.array(range(1,self.N+1))/ (self.N+1)
  
class TestSolver(unittest.TestCase):
  """
    The way this works it tests whichever nonlinear solver was built into 
    the system.
  """
  def setUp(self):
    self.N = 15
    self.model = solvers.TestRosenbrock(self.N)
    self.soln = np.ones((self.N,))

  def test_solve(self):
    dummy = solvers.TrialState()
    x = solvers.solve(self.model, dummy)
    self.assertTrue(np.allclose(self.soln, x))
