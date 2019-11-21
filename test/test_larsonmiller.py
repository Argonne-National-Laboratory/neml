import sys
sys.path.append('..')

from neml import larsonmiller, interpolate
import unittest

from common import *

import numpy as np

class TestLarsonMiller(unittest.TestCase):
  def setUp(self):
    self.fn = interpolate.PolynomialInterpolate([-6.653e-9,2.952e-4,-6.197e-1])
    self.C = 32.06

    self.lmr = larsonmiller.LarsonMillerRelation(self.fn, self.C)

    self.T = 550 + 273.15

  def test_sR(self):
    tR = 50000.0
    lmp = self.T*(self.C + np.log10(tR))
    direct = 10.0**(self.fn.value(lmp))
    test = self.lmr.sR(tR, self.T)
    self.assertTrue(np.isclose(direct, test))

  def test_tR(self):
    tR = 50000.0
    lmp = self.T*(self.C + np.log10(tR))
    sR = 10.0**(self.fn.value(lmp))
    tR2 = self.lmr.tR(sR, self.T)
    self.assertTrue(np.isclose(tR2, tR))

  def test_dtR(self):
    sR = 200.0
    num = differentiate(lambda s: self.lmr.tR(s, self.T), sR)
    act = self.lmr.dtR_ds(sR, self.T)
    self.assertTrue(np.isclose(num, act, rtol = 1e-4))
