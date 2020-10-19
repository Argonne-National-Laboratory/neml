import sys
sys.path.append('..')

import unittest

import numpy as np

from neml import solvers

class TestConvergence(unittest.TestCase):
  def setUp(self):
    self.A = 1.0
    self.n = 4.0
    self.b = -16.0
    self.x0 = 4.0
    self.model = solvers.TestPower(self.A, self.n, self.b, self.x0)
    self.ts = solvers.TrialState()

  def test_correct(self):
    x = solvers.solve(self.model, self.ts, rtol = 1.0e-20, atol = 1.0e-12,
        linesearch = False, verbose = False)
    self.assertAlmostEqual(x[0], (-self.b/self.A)**(1.0/self.x0))

  def test_rtol(self):
    R0, J0 = self.model.RJ(np.array([self.x0]), self.ts)

    x = solvers.solve(self.model, self.ts, rtol = 1.0e-4, atol = 1.0e-40,
        linesearch = False, verbose = False)
    Rf, Rj = self.model.RJ(x, self.ts)
    self.assertTrue(Rf[0]/R0[0] < 1.0e-4)

  def test_atol(self):
    R0, J0 = self.model.RJ(np.array([self.x0]), self.ts)

    x = solvers.solve(self.model, self.ts, rtol = 1.0e-40, atol = 1.0e-4,
        linesearch = False, verbose = False)
    Rf, Rj = self.model.RJ(x, self.ts)
    self.assertTrue(Rf[0] < 1.0e-4)
