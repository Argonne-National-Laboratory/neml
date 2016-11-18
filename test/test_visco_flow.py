import sys
sys.path.append('..')

from neml import visco_flow, surfaces, hardening
from common import *

import unittest
import numpy as np
import numpy.linalg as la
import numpy.random as ra

class CommonGFlow(object):
  """
    Common to our Perzyna g functions.
  """
  def test_dg(self):
    s = self.gen_f()

    exact = self.model.dg(s)

    dfn = lambda x: self.model.g(x)
    num = differentiate(dfn, s)

    self.assertTrue(np.isclose(exact, num))

class TestGPowerLaw(unittest.TestCase, CommonGFlow):
  def setUp(self):
    self.n = 11.0
    self.model = visco_flow.GPowerLaw(self.n)

  def gen_f(self):
    return 100*(1-2.0*ra.random((1,))[0])

  def test_properties(self):
    self.assertTrue(np.isclose(self.n, self.model.n))

  def test_g(self):
    s = self.gen_f()
    if (s > 0.0):
      v = s**self.n
    else:
      v = 0.0
    self.assertTrue(np.isclose(self.model.g(s), v))


class CommonFlowRule(object):
  """
    Tests common to all flow rules.
  """
  def gen_stress(self):
    s = np.array([200,0,100,0,-50,0])
    return s

  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_dy_ds(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.y(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dy_ds(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3, atol = 1.0e-3))

  def test_dy_da(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.y(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dy_da(stress, hist, self.T)
  
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_ds(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.g(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dg_ds(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_da(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.g(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dg_da(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_ds_time(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.g_time(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dg_ds_time(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_da_time(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.g_time(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dg_da_time(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_ds_temp(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.g_temp(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dg_ds_temp(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dg_da_temp(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.g_temp(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dg_da_temp(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_ds(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.h(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dh_ds(stress, hist, self.T)
    
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_da(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.h(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dh_da(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_ds_time(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.h_time(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dh_ds_time(stress, hist, self.T)
    
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_da_time(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.h_time(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dh_da_time(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_ds_temp(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda s: self.model.h_temp(s, hist, self.T)
    num = differentiate(dfn, stress)
    exact = self.model.dh_ds_temp(stress, hist, self.T)
    
    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))

  def test_dh_da_temp(self):
    stress = self.gen_stress()
    hist = self.gen_hist()

    dfn = lambda a: self.model.h_temp(stress, a, self.T)
    num = differentiate(dfn, hist)
    exact = self.model.dh_da_temp(stress, hist, self.T)

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3))


class TestPerzynaIsoJ2Voce(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.s0 = 180.0
    self.R = 150.0
    self.d = 10.0

    self.n = 5.0
    self.eta = 20.0

    self.surface = surfaces.IsoJ2()
    self.hrule = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)
    self.g = visco_flow.GPowerLaw(self.n)

    self.model = visco_flow.PerzynaFlowRule(self.surface, self.hrule, self.g, self.eta)

    self.hist0 = np.zeros((1,))
    self.T = 300.0

  def gen_hist(self):
    return np.array([0.01])

  def test_properties(self):
    self.assertTrue(np.isclose(self.eta, self.model.eta))


class TestPerzynaIsoJ2Linear(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.s0 = 100.0
    self.Kp = 1500.0

    self.n = 2.0
    self.eta = 100.0

    self.surface = surfaces.IsoJ2()
    self.hrule = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    self.g = visco_flow.GPowerLaw(self.n)

    self.model = visco_flow.PerzynaFlowRule(self.surface, self.hrule, self.g, self.eta)

    self.hist0 = np.zeros((1,))
    self.T = 300.0

  def gen_hist(self):
    return np.array([0.01])

  def test_properties(self):
    self.assertTrue(np.isclose(self.eta, self.model.eta))

class TestChebocheModel(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    n = 20.0
    eta = 108.0
    sY = 89.0

    Q = 165.0
    b = 12.0
    
    self.m = 3

    C1 = 80.0e3
    C2 = 14.02e3
    C3 = 3.333e3

    y1 = 0.9e3
    y2 = 1.5e3
    y3 = 1.0

    surface = surfaces.IsoKinJ2()
    iso = hardening.VoceIsotropicHardeningRule(sY, Q, b)
    cs = np.array([C1, C2, C3])
    gs = np.array([y1, y2, y3])
    hmodel = hardening.Chaboche(iso, cs, gs)

    fluidity = visco_flow.ConstantFluidity(eta)

    self.hist0 = np.zeros((19,))
    self.T = 300.0

    self.model = visco_flow.ChabocheFlowRule(
        surface, hmodel, fluidity, n)

  def gen_hist(self):
    h = np.array([0.05,20,-30,40.0,5.0,2.0,40.0,-10,-20,10,50,30,10,-50,60,70,15,-15,10])
    h[1:7] = make_dev(h[1:7])
    h[7:13] = make_dev(h[7:13])
    h[13:19] = make_dev(h[13:19])
    return h

class TestChebocheFlow(unittest.TestCase):
  def setUp(self):
    self.n = 20.0
    self.eta = 108.0
    self.sY = 89.0

    self.Q = 165.0
    self.b = 12.0
    
    self.m = 3

    C1 = 80.0e3
    C2 = 14.02e3
    C3 = 3.333e3

    y1 = 0.9e3
    y2 = 1.5e3
    y3 = 1.0

    surface = surfaces.IsoKinJ2()
    self.iso = hardening.VoceIsotropicHardeningRule(self.sY, self.Q, self.b)
    self.cs = np.array([C1, C2, C3])
    self.gs = np.array([y1, y2, y3])
    self.hmodel = hardening.Chaboche(self.iso, self.cs, self.gs)

    self.fluidity = visco_flow.ConstantFluidity(self.eta)

    self.hist0 = np.zeros((19,))
    self.T = 300.0

    self.model = visco_flow.ChabocheFlowRule(
        surface, self.hmodel, self.fluidity, self.n)

  def gen_stress(self):
    s = (2.0 * ra.random((6,)) - 1.0) * 400
    return s

  def gen_hist(self):
    h = ra.random((self.m*6+1,))
    h[1:self.m*6+1] = (2.0 * h[1:self.m*6+1] - 1.0) * 200
    for i in range(self.m):
      h[1+i*6:1+(i+1)*6] = make_dev(h[1+i*6:1+(i+1)*6])
    return h

  def test_y(self):
    s = self.gen_stress()
    h = self.gen_hist()
    y = self.model.y(s, h, self.T)

    sdev = make_dev(s)
    X = sum(h[1+i*6:1+(i+1)*6] for i in range(self.m))
    ih = self.iso.q(h[0:1], self.T)[0]

    tv = np.sqrt(3.0/2.0) * la.norm(sdev + X) + ih

    if tv > 0.0:
      pd = (tv / self.eta)**self.n
    else:
      pd = 0.0

    gd = pd / np.sqrt(2.0/3.0)

    self.assertTrue(np.isclose(gd, y))

  def test_ep(self):
    s = self.gen_stress()
    h = self.gen_hist()
    y = self.model.y(s, h, self.T)

    ep_model = self.model.g(s, h, self.T) * y

    sdev = make_dev(s)
    X = sum(h[1+i*6:1+(i+1)*6] for i in range(self.m))
    ih = self.iso.q(h[0:1], self.T)[0]

    tv = np.sqrt(3.0/2.0) * la.norm(sdev + X) + ih

    if tv > 0.0:
      pd = (tv / self.eta)**self.n
    else:
      pd = 0.0

    ep_calc = 3.0/2.0 * pd * (sdev + X) / (np.sqrt(3.0/2.0) * la.norm(sdev+X))

    self.assertTrue(np.allclose(ep_calc, ep_model))

  def test_rdot(self):
    s = self.gen_stress()
    h = self.gen_hist()
    y = self.model.y(s, h, self.T)
    adot = self.model.h(s,h,self.T)[0] * y
    dr = -self.b * np.exp(-self.b * h[0]) * self.Q * adot

    sdev = make_dev(s)
    X = sum(h[1+i*6:1+(i+1)*6] for i in range(self.m))
    ih = self.iso.q(h[0:1], self.T)[0]

    tv = np.sqrt(3.0/2.0) * la.norm(sdev + X) + ih

    if tv > 0.0:
      pd = (tv / self.eta)**self.n
    else:
      pd = 0.0

    q = self.hmodel.q(h, self.T)

    Rtot = -self.sY - self.Q * (1.0 - np.exp(-self.b*h[0]))

    self.assertTrue(np.isclose(Rtot, q[0]))

    R_act = -(Rtot + self.sY)

    dr_direct = -self.b*(self.Q - R_act) * pd

    self.assertTrue(np.isclose(dr_direct, dr))

  def test_Xdot(self):
    s = self.gen_stress()
    h = self.gen_hist()
    y = self.model.y(s, h, self.T)
    dX = self.model.h(s, h, self.T)[1:] * y

    sdev = make_dev(s)
    X = sum(h[1+i*6:1+(i+1)*6] for i in range(self.m))
    ih = self.iso.q(h[0:1], self.T)[0]

    tv = np.sqrt(3.0/2.0) * la.norm(sdev + X) + ih

    if tv > 0.0:
      pd = (tv / self.eta)**self.n
    else:
      pd = 0.0

    ep_calc = 3.0/2.0 * pd * (sdev + X) / (np.sqrt(3.0/2.0) * la.norm(sdev+X))

    dX_calc = np.zeros((6*self.m,))

    for i in range(self.m):
      dX_calc[i*6:(i+1)*6] = -2.0/3.0 * self.cs[i] * ep_calc - self.gs[i] * h[1+i*6:1+(i+1)*6] * pd

    self.assertTrue(np.allclose(dX, dX_calc))


class CommonFluidity(object):
  """
    Common fluidity model tests
  """
  def gen_hist(self):
    return ra.random((1,))[0]

  def test_deta(self):
    a = self.gen_hist()
    fm = self.model.deta(a)
    dfn = lambda a: self.model.eta(a)
    fn = differentiate(dfn, a)
    self.assertTrue(np.isclose(fm, fn))

class TestConstantFluidity(unittest.TestCase, CommonFluidity):
  def setUp(self):
    self.eta = 200.0
    self.model = visco_flow.ConstantFluidity(self.eta)

  def test_eta(self):
    a = self.gen_hist()
    self.assertTrue(np.isclose(self.eta, self.model.eta(a)))

class TestChabocheJ2Voce(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.n = 12.0
    self.K = 150.0
    self.k = 6.0
    self.C = 24800.0
    self.g0 = 300.0
    self.Q = 86 - self.k
    self.gs = 300.0
    self.b = 10.0
    self.beta = 0.0

    self.surface = surfaces.IsoKinJ2()

    self.iso = hardening.VoceIsotropicHardeningRule(self.k, self.Q, self.b)
    cs = [self.C]
    gs = [hardening.SatGamma(self.gs, self.g0, self.beta)]
    self.hardening = hardening.Chaboche(self.iso, cs, gs)

    self.fluidity = visco_flow.ConstantFluidity(self.K)

    self.model = visco_flow.ChabocheFlowRule(
        self.surface, self.hardening, self.fluidity, self.n) 

    self.hist0 = np.zeros((7,))
    self.T = 300.0

  def gen_hist(self):
    bs = make_dev([10.0, -50.0, 25.0, 30.0, -10.0, 30.0])
    return np.array([0.01] + list(bs))

