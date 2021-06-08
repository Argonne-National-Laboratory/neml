import sys
sys.path.append('..')

from neml import visco_flow, surfaces, hardening
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonGFlow(object):
  """
    Common to our Perzyna g functions.
  """
  def test_dg(self):
    s = self.gen_f()

    exact = self.model.dg(s, self.T)

    dfn = lambda x: self.model.g(x, self.T)
    num = differentiate(dfn, s)

    self.assertTrue(np.isclose(exact, num))

class TestGPowerLaw(unittest.TestCase, CommonGFlow):
  def setUp(self):
    self.n = 11.0
    self.eta = 200.0
    self.model = visco_flow.GPowerLaw(self.n, self.eta)
    self.T = 300.0

  def gen_f(self):
    return 100*0.25

  def test_properties(self):
    self.assertTrue(np.isclose(self.n, self.model.n(self.T)))
    self.assertTrue(np.isclose(self.eta, self.model.eta(self.T)))

  def test_g(self):
    s = self.gen_f()
    if (s > 0.0):
      v = (s/self.eta)**self.n
    else:
      v = 0.0
    self.assertTrue(np.isclose(self.model.g(s, self.T), v))

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

    self.assertTrue(np.allclose(num, exact, rtol = 1.0e-3, atol = 1.0e-6))

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
    self.g = visco_flow.GPowerLaw(self.n, self.eta)

    self.model = visco_flow.PerzynaFlowRule(self.surface, self.hrule, self.g)

    self.hist0 = np.zeros((1,))
    self.T = 300.0

  def gen_hist(self):
    return np.array([0.01])

class TestPerzynaIsoJ2Linear(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.s0 = 100.0
    self.Kp = 1500.0

    self.n = 2.0
    self.eta = 100.0

    self.surface = surfaces.IsoJ2()
    self.hrule = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    self.g = visco_flow.GPowerLaw(self.n, self.eta)

    self.model = visco_flow.PerzynaFlowRule(self.surface, self.hrule, self.g)

    self.hist0 = np.zeros((1,))
    self.T = 300.0

  def gen_hist(self):
    return np.array([0.01])

class TestChabocheModel(unittest.TestCase, CommonFlowRule):
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
    cs = [C1, C2, C3]
    gs = [y1, y2, y3]
    gmodels = [hardening.ConstantGamma(g) for g in gs]
    A = [0.0, 0.0, 0.0]
    ae = [1.0, 1.0, 1.0]

    hmodel = hardening.Chaboche(iso, cs, gmodels, A, ae)

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

class TestChabocheFlow(unittest.TestCase):
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
    self.cs = [C1, C2, C3]
    self.gs = [y1, y2, y3]
    self.As = [0.0, 0.0, 0.0]
    self.ns = [1.0, 1.0, 1.0]
    self.gmodels = [hardening.ConstantGamma(g) for g in self.gs]
    self.hmodel = hardening.Chaboche(self.iso, self.cs, self.gmodels, 
        self.As, self.ns)

    self.fluidity = visco_flow.ConstantFluidity(self.eta)

    self.hist0 = np.zeros((19,))
    self.T = 300.0

    self.model = visco_flow.ChabocheFlowRule(
        surface, self.hmodel, self.fluidity, self.n)

  def gen_stress(self):
    return np.array([300.0,-200.0,100.0,50.0,150.0,300.0])

  def gen_hist(self):
    h = np.array(range(1,self.m*6+2)) / (self.m*7)
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


class TestChabocheFlowWithPrefactor(unittest.TestCase):
  def setUp(self):
    self.n = 20.0
    self.eta = 108.0
    self.sY = 89.0
    self.prefactor = np.asarray([2.])
    self.prefactor = 2.

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
    self.cs = [C1, C2, C3]
    self.gs = [y1, y2, y3]
    self.As = [0.0, 0.0, 0.0]
    self.ns = [1.0, 1.0, 1.0]
    self.gmodels = [hardening.ConstantGamma(g) for g in self.gs]
    self.hmodel = hardening.Chaboche(self.iso, self.cs, self.gmodels,
        self.As, self.ns)

    self.fluidity = visco_flow.ConstantFluidity(self.eta)

    self.hist0 = np.zeros((19,))
    self.T = 300.0

    self.model = visco_flow.ChabocheFlowRule(
        surface, self.hmodel, self.fluidity, self.n, prefactor = self.prefactor)

  def gen_stress(self):
    return np.array([300.0,-200.0,100.0,50.0,150.0,300.0])

  def gen_hist(self):
    h = np.array(range(1,self.m*6+2)) / (self.m*7)
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
      pd = self.prefactor*(tv / self.eta)**self.n
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
      pd = self.prefactor*(tv / self.eta)**self.n
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
      pd = self.prefactor*(tv / self.eta)**self.n
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
      pd = self.prefactor*(tv / self.eta)**self.n
    else:
      pd = 0.0

    ep_calc = 3.0/2.0 * pd * (sdev + X) / (np.sqrt(3.0/2.0) * la.norm(sdev+X))

    dX_calc = np.zeros((6*self.m,))

    for i in range(self.m):
      dX_calc[i*6:(i+1)*6] = -2.0/3.0 * self.cs[i] * ep_calc - self.gs[i] * h[1+i*6:1+(i+1)*6] * pd

    self.assertTrue(np.allclose(dX, dX_calc))

class TestChabocheModelWithFactor(unittest.TestCase, CommonFlowRule):
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
    cs = [C1, C2, C3]
    gs = [y1, y2, y3]
    gmodels = [hardening.ConstantGamma(g) for g in gs]
    A = [0.0, 0.0, 0.0]
    ae = [1.0, 1.0, 1.0]

    self.prefactor = 2.0

    hmodel = hardening.Chaboche(iso, cs, gmodels, A, ae)

    fluidity = visco_flow.ConstantFluidity(eta)

    self.hist0 = np.zeros((19,))
    self.T = 300.0

    self.model = visco_flow.ChabocheFlowRule(
        surface, hmodel, fluidity, n, prefactor = self.prefactor)

  def gen_hist(self):
    h = np.array([0.05,20,-30,40.0,5.0,2.0,40.0,-10,-20,10,50,30,10,-50,60,70,15,-15,10])
    h[1:7] = make_dev(h[1:7])
    h[7:13] = make_dev(h[7:13])
    h[13:19] = make_dev(h[13:19])
    return h

class CommonFluidity(object):
  """
    Common fluidity model tests
  """
  def gen_hist(self):
    return 0.15

  def test_deta(self):
    a = self.gen_hist()
    fm = self.model.deta(a, self.T)
    dfn = lambda a: self.model.eta(a, self.T)
    fn = differentiate(dfn, a)
    self.assertTrue(np.isclose(fm, fn))

class TestConstantFluidity(unittest.TestCase, CommonFluidity):
  def setUp(self):
    self.eta = 200.0
    self.model = visco_flow.ConstantFluidity(self.eta)
    self.T = 300.0

  def test_eta(self):
    a = self.gen_hist()
    self.assertTrue(np.isclose(self.eta, self.model.eta(a, self.T)))

class TestSaturatingFluidity(unittest.TestCase, CommonFluidity):
  def setUp(self):
    self.K0 = 116.0
    self.A = 238.9
    self.b = 9.80
    self.model = visco_flow.SaturatingFluidity(self.K0, self.A, self.b)
    self.T = 300.0

  def test_eta(self):
    a = self.gen_hist()
    self.assertTrue(np.isclose(self.K0 + self.A * (1.0 - np.exp(-self.b * a)), self.model.eta(a, self.T)))

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
    As = [0.0]
    ns = [1.0]
    self.hardening = hardening.Chaboche(self.iso, cs, gs, As, ns)

    self.fluidity = visco_flow.ConstantFluidity(self.K)

    self.model = visco_flow.ChabocheFlowRule(
        self.surface, self.hardening, self.fluidity, self.n) 

    self.hist0 = np.zeros((7,))
    self.T = 300.0

  def gen_hist(self):
    bs = make_dev([10.0, -50.0, 25.0, 30.0, -10.0, 30.0])
    return np.array([0.01] + list(bs))

class TestYaguchiGr91Flow(unittest.TestCase, CommonFlowRule):
  def setUp(self):
    self.model = visco_flow.YaguchiGr91FlowRule()

    self.pTrange = np.linspace(473.0, 873.0)

    self.hist0 = np.zeros((14,))
    self.T = 500.0

  def gen_stress(self):
    s = np.array([200,75,100,50,-50,100])*0.1
    return s
  
  @staticmethod
  def J2(x):
    return np.sqrt(3.0/2.0 * np.dot(x,x))

  def test_y(self):
    s = self.gen_stress()
    h = self.gen_hist()

    y_model = self.model.y(s, h, self.T)

    X = h[:6] + h[6:12]
    Xdev = make_dev(X)
    Sdev = make_dev(s)

    t1 = (self.J2(Xdev - Sdev) - h[13])/self.model.D(self.T)
    if t1 > 0.0:
      y_calc = t1**self.model.n(self.T)
    else:
      y_calc = 0.0

    self.assertTrue(np.isclose(y_calc, y_model))

  def test_g(self):
    s = self.gen_stress()
    h = self.gen_hist()

    g_model = self.model.g(s, h, self.T)

    X = h[:6] + h[6:12]
    Xdev = make_dev(X)
    Sdev = make_dev(s)

    g_calc = 3.0/2.0 * (Sdev - Xdev) / self.J2(Sdev - Xdev)

    self.assertTrue(np.allclose(g_calc, g_model))

  def test_h(self):
    s = self.gen_stress()
    a = self.gen_hist()

    h_model = self.model.h(s, a, self.T)
    h_calc = np.zeros(h_model.shape)
    
    C1 = self.model.C1(self.T)
    C2 = self.model.C2(self.T)

    a1 = self.model.a10(self.T) - a[12]
    a2 = self.model.a2(self.T)

    X = a[:6] + a[6:12]
    Xdev = make_dev(X)
    Sdev = make_dev(s)

    n = 3.0/2.0 * (Sdev - Xdev) / self.J2(Sdev - Xdev)

    h_calc[:6] = C1 * (2.0/3.0 * a1 * n - a[:6])

    h_calc[6:12] = C2 * (2.0/3.0 * a2 * n - a[6:12])

    h_calc[12] = self.model.d(self.T) * (self.model.q(self.T) - a[12])

    t1 = (self.J2(Xdev - Sdev) - a[13])/self.model.D(self.T)
    if t1 > 0.0:
      y_calc = t1**self.model.n(self.T)
    else:
      y_calc = 0.0

    sas = self.model.A(self.T) + self.model.B(self.T) * np.log10(y_calc)
    if sas < 0.0:
      sas = 0.0

    if (sas - a[13]) >= 0:
      b = self.model.bh(self.T)
    else:
      b = self.model.br(self.T)

    h_calc[13] = b * (sas - a[13])

    self.assertTrue(np.allclose(h_calc, h_model))

  def gen_hist(self):
    X1 = make_dev([15.0,-20.0,-30.0,50.0,-10.0,5])
    X2 = make_dev([-25,10,15,30,-15,20])
    return np.array(list(X1) + list(X2) + [50.0] + [5.0]) / 10.0

  def eval_prop(self, prop):
    return np.array([prop(T) for T in self.pTrange])

  def test_D(self):
    Ds_model = self.eval_prop(self.model.D)
    Ds_direct = (-1.119e-9 + 5.145e-11 * self.pTrange - 5.450e-14 * self.pTrange**2.0)**(-1.0/2.85)
    self.assertTrue(np.allclose(Ds_model, Ds_direct))
  
  def test_n(self):
    ns_model = self.eval_prop(self.model.n)
    ns_direct = 2.850 * np.ones(self.pTrange.shape)
    self.assertTrue(np.allclose(ns_model, ns_direct))

  def test_a10(self):
    a10s_model = self.eval_prop(self.model.a10)
    a10s_direct = 2.082e3 - 8.110*self.pTrange + 1.321e-2*self.pTrange**2.0 - 7.278e-6 * self.pTrange**3.0
    self.assertTrue(np.allclose(a10s_model, a10s_direct))

  def test_C2(self):
    C2s_model = self.eval_prop(self.model.C2)
    C2s_direct = np.piecewise(self.pTrange, [self.pTrange < 798.0, self.pTrange >= 798.0], [lambda x: 2.000e2, lambda x: -2.992e3 + 4.0 * x])
    self.assertTrue(np.allclose(C2s_model, C2s_direct))

  def test_a2(self):
    a2s_model = self.eval_prop(self.model.a2)
    a2s_direct = np.piecewise(self.pTrange, [self.pTrange < 773.0, self.pTrange >= 773.0], [lambda x: 1.100e2, lambda x: 3.373e3 - 7.622*x + 4.400e-3*x**2.0])
    self.assertTrue(np.allclose(a2s_model, a2s_direct))

  def test_g1(self):
    g1s_model = self.eval_prop(self.model.g1)
    g1s_direct = np.piecewise(self.pTrange, [self.pTrange < 773.0, self.pTrange >= 773.0], [lambda x: 5.064e-137 * np.exp(3.545e-1 * x), lambda x: 5.336e-33 * np.exp(4.470e-2 * x)])
    self.assertTrue(np.allclose(g1s_model, g1s_direct))

  def test_g2(self):
    g2s_model = self.eval_prop(self.model.g2)
    g2s_direct = np.piecewise(self.pTrange, [self.pTrange < 773.0, self.pTrange >= 773.0], [lambda x: 8.572e-108 * np.exp(2.771e-1*x), lambda x: 2.817e-11 - 7.538e-14 * x + 5.039e-17 * x**2.0])
    self.assertTrue(np.allclose(g2s_model, g2s_direct))
  
  def test_m(self):
    ms_model = self.eval_prop(self.model.m)
    ms_direct = np.piecewise(self.pTrange, [self.pTrange < 673.0, np.logical_and(673.0 <= self.pTrange, self.pTrange < 773.0), self.pTrange >= 773.0], [lambda x: 1.200e1, lambda x: 5.766e2*np.exp(-5.754e-3*x), lambda x: 6.750])
    self.assertTrue(np.allclose(ms_model, ms_direct))
  
  def test_br(self):
    brs_model = self.eval_prop(self.model.br)
    brs_direct = 1.000e3 * np.ones(self.pTrange.shape)
    self.assertTrue(np.allclose(brs_model, brs_direct))

  def test_bh(self):
    bhs_model = self.eval_prop(self.model.bh)
    bhs_direct = 1.065e3 - 1.404 * self.pTrange + 4.417e-4 * self.pTrange**2.0
    self.assertTrue(np.allclose(bhs_model, bhs_direct))

  def test_A(self):
    As_model = self.eval_prop(self.model.A)
    As_direct = np.piecewise(self.pTrange, [self.pTrange < 673.0, np.logical_and(673.0 <= self.pTrange, self.pTrange < 823.0), self.pTrange >= 823.0], [lambda x: -1.841e2 + 2.075e-1 * x, lambda x: -4.799e2 + 1.262*x - 9.133e-4 * x**2.0, lambda x: -6.0e1])
    self.assertTrue(np.allclose(As_model, As_direct))

  def test_B(self):
    Bs_model = self.eval_prop(self.model.B)
    Bs_direct = np.piecewise(self.pTrange, [self.pTrange < 773.0, self.pTrange >= 773.0], [lambda x: -6.340e1 + 6.850e-2*x, lambda x: -1.730e1])
    self.assertTrue(np.allclose(Bs_model, Bs_direct))

  def test_d(self):
    ds_model = self.eval_prop(self.model.d)
    ds_direct = 3.208 - 1.010e-2 * self.pTrange + 1.012e-5 * self.pTrange**2.0
    self.assertTrue(np.allclose(ds_model, ds_direct))

  def test_q(self):
    qs_model = self.eval_prop(self.model.q)
    qs_direct = np.piecewise(self.pTrange, [self.pTrange < 673.0, np.logical_and(673.0 <= self.pTrange, self.pTrange < 823.0), self.pTrange >= 823.0], [lambda x: 9.50e1, lambda x: -2.467e2 + 7.320e-1*x - 3.333e-4*x**2.0, lambda x: 1.300e2])
    self.assertTrue(np.allclose(qs_model, qs_direct))

  def test_C1(self):
    C1s_model = self.eval_prop(self.model.C1)
    C1s_direct = np.piecewise(self.pTrange, [self.pTrange < 673.0, np.logical_and(673.0 <= self.pTrange, self.pTrange < 773.0), self.pTrange >= 773.0], [lambda x: 1.50e3, lambda x: -2.879e4 + 45.0*x, lambda x: 6.000e3])
    self.assertTrue(np.allclose(C1s_model, C1s_direct))

