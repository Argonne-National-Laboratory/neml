import sys
sys.path.append('..')

from neml import general_flow, visco_flow, elasticity, surfaces, hardening
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonGeneralFlow(object):
  def test_ds_ds(self):
    t_np1 = self.gen_t()
    e_np1 = self.gen_e()
    e_dot = self.gen_edot(e_np1, t_np1)
    T_np1 = self.gen_T()
    T_dot = self.gen_Tdot(T_np1, t_np1)
    s_np1 = self.gen_stress()
    h_np1 = self.gen_hist()

    dfn = lambda x: self.model.s(x, h_np1, e_dot, T_np1, T_dot)
    num = differentiate(dfn, s_np1)
    should = self.model.ds_ds(s_np1, h_np1, e_dot, T_np1, T_dot)
    
    self.assertTrue(np.allclose(num, should))

  def test_ds_da(self):
    t_np1 = self.gen_t()
    e_np1 = self.gen_e()
    e_dot = self.gen_edot(e_np1, t_np1)
    T_np1 = self.gen_T()
    T_dot = self.gen_Tdot(T_np1, t_np1)
    s_np1 = self.gen_stress()
    h_np1 = self.gen_hist()

    dfn = lambda x: self.model.s(s_np1, x, e_dot, T_np1, T_dot)
    num = differentiate(dfn, h_np1)
    should = self.model.ds_da(s_np1, h_np1, e_dot, T_np1, T_dot)
    
    self.assertTrue(np.allclose(num, should, rtol = 1.0e-3))

  def test_ds_de(self):
    t_np1 = self.gen_t()
    e_np1 = self.gen_e()
    e_dot = self.gen_edot(e_np1, t_np1)
    T_np1 = self.gen_T()
    T_dot = self.gen_Tdot(T_np1, t_np1)
    s_np1 = self.gen_stress()
    h_np1 = self.gen_hist()

    dfn = lambda x: self.model.s(s_np1, h_np1, x, T_np1, T_dot)
    num = differentiate(dfn, e_dot, eps = self.eps)
    should = self.model.ds_de(s_np1, h_np1, e_dot, T_np1, T_dot)

    self.assertTrue(np.allclose(num, should, rtol = 1.0e-3))

  def test_da_ds(self):
    t_np1 = self.gen_t()
    e_np1 = self.gen_e()
    e_dot = self.gen_edot(e_np1, t_np1)
    T_np1 = self.gen_T()
    T_dot = self.gen_Tdot(T_np1, t_np1)
    s_np1 = self.gen_stress()
    h_np1 = self.gen_hist()

    dfn = lambda x: self.model.a(x, h_np1, e_dot, T_np1, T_dot)
    num = differentiate(dfn, s_np1)
    should = self.model.da_ds(s_np1, h_np1, e_dot, T_np1, T_dot)
    
    self.assertTrue(np.allclose(num, should, rtol = 1.0e-3))

  def test_da_da(self):
    t_np1 = self.gen_t()
    e_np1 = self.gen_e()
    e_dot = self.gen_edot(e_np1, t_np1)
    T_np1 = self.gen_T()
    T_dot = self.gen_Tdot(T_np1, t_np1)
    s_np1 = self.gen_stress()
    h_np1 = self.gen_hist()

    dfn = lambda x: self.model.a(s_np1, x, e_dot, T_np1, T_dot)
    num = differentiate(dfn, h_np1)
    should = self.model.da_da(s_np1, h_np1, e_dot, T_np1, T_dot)

    self.assertTrue(np.allclose(num, should, rtol = 1.0e-3))

  def test_da_de(self):
    t_np1 = self.gen_t()
    e_np1 = self.gen_e()
    e_dot = self.gen_edot(e_np1, t_np1)
    T_np1 = self.gen_T()
    T_dot = self.gen_Tdot(T_np1, t_np1)
    s_np1 = self.gen_stress()
    h_np1 = self.gen_hist()

    dfn = lambda x: self.model.a(s_np1, h_np1, x, T_np1, T_dot)
    num = differentiate(dfn, e_dot, eps = self.eps)
    should = self.model.da_de(s_np1, h_np1, e_dot, T_np1, T_dot)

    self.assertTrue(np.allclose(num, should, rtol = 1.0e-3))


class CommonTVPFlow(object):
  def test_history(self):
    self.assertEqual(len(self.h_n), self.model.nhist)
    self.assertTrue(np.allclose(self.h_n, self.model.init_hist()))

  def test_srate(self):
    t_np1 = self.gen_t()
    e_np1 = self.gen_e()
    e_dot = self.gen_edot(e_np1, t_np1)
    T_np1 = self.gen_T()
    T_dot = self.gen_Tdot(T_np1, t_np1)
    s_np1 = self.gen_stress()
    h_np1 = self.gen_hist()

    should = np.dot(self.emodel.C(T_np1), 
        e_dot - self.vmodel.g(s_np1, h_np1, T_np1) * 
        self.vmodel.y(s_np1, h_np1, T_np1) - 
        self.vmodel.g_temp(s_np1, h_np1, T_np1) * T_dot -
        self.vmodel.g_time(s_np1, h_np1, T_np1))

    actual = self.model.s(s_np1, h_np1, e_dot, T_np1, T_dot)

    self.assertTrue(np.allclose(should, actual))

  def test_hrate(self):
    t_np1 = self.gen_t()
    e_np1 = self.gen_e()
    e_dot = self.gen_edot(e_np1, t_np1)
    T_np1 = self.gen_T()
    T_dot = self.gen_Tdot(T_np1, t_np1)
    s_np1 = self.gen_stress()
    h_np1 = self.gen_hist()

    should = (self.vmodel.h(s_np1, h_np1, T_np1) * 
      self.vmodel.y(s_np1, h_np1, T_np1) + 
      self.vmodel.h_temp(s_np1, h_np1, T_np1) * T_dot +
      self.vmodel.h_time(s_np1, h_np1, T_np1))

    actual = self.model.a(s_np1, h_np1, e_dot, T_np1, T_dot)

    self.assertTrue(np.allclose(should, actual))


class TestTVPCheboche(unittest.TestCase, CommonGeneralFlow, CommonTVPFlow):
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
    As = [0.0, 0.0, 0.0]
    ns = [1.0, 1.0, 1.0]
    gmodels = [hardening.ConstantGamma(g) for g in gs]
    hmodel = hardening.Chaboche(iso, cs, gmodels, As, ns)

    fluidity = visco_flow.ConstantFluidity(eta)

    self.vmodel = visco_flow.ChabocheFlowRule(
        surface, hmodel, fluidity, n)

    E = 92000.0
    nu = 0.3

    mu = E/(2*(1+nu))
    K = E/(3*(1-2*nu))

    self.emodel = elasticity.IsotropicLinearElasticModel(mu, "shear",
        K, "bulk")

    self.model = general_flow.TVPFlowRule(self.emodel, self.vmodel)

    self.T_n = 300.0
    self.e_n = np.zeros((6,))
    self.t_n = 0.0
    self.h_n = np.zeros((1+6*self.m,))

    self.eps = 1.0e-3

  def gen_hist(self):
    h = np.array([0.1,50,-30,40,60,-80,-30,-100,-15,-30,-50,100,50,60,-60,50,
      30,90,40])
    for i in range(self.m):
      h[1+i*6:1+(i+1)*6] = make_dev(h[1+i*6:1+(i+1)*6])
    
    return h

  def gen_stress(self):
    return np.array([200.0,-200.0,100.0,50.0,25.0,-50.0])

  def gen_e(self):
    return np.array([0.025,-0.01,-0.02,0.01,0.02,-0.03])

  def gen_T(self):
    return 350

  def gen_t(self):
    return 1.25

  def gen_dt(self, t):
    return t - self.t_n

  def gen_edot(self, e, t):
    return (e - self.e_n) / self.gen_dt(t)

  def gen_Tdot(self, T, t):
    return (T - self.T_n) / self.gen_dt(t)


class TestTVPYaguchi(unittest.TestCase, CommonGeneralFlow, CommonTVPFlow):
  def setUp(self):
    self.vmodel = visco_flow.YaguchiGr91FlowRule()

    E = 92000.0
    nu = 0.3

    mu = E/(2*(1+nu))
    K = E/(3*(1-2*nu))

    self.emodel = elasticity.IsotropicLinearElasticModel(mu, "shear",
        K, "bulk")

    self.model = general_flow.TVPFlowRule(self.emodel, self.vmodel)

    self.T_n = 300.0
    self.e_n = np.zeros((6,))
    self.t_n = 0.0
    self.h_n = np.zeros((2+12,))

    self.eps = 1.0e-3

  def gen_hist(self):
    h = np.array([50,-30,40,60,-80,-30,-100,-15,-30,-50,100,50,5,2.5])
    for i in range(2):
      h[1+i*6:1+(i+1)*6] = make_dev(h[1+i*6:1+(i+1)*6])
    
    return h

  def gen_stress(self):
    return np.array([200.0,-200.0,100.0,50.0,25.0,-50.0])

  def gen_e(self):
    return np.array([0.025,-0.01,-0.02,0.01,0.02,-0.03])

  def gen_T(self):
    return 350

  def gen_t(self):
    return 1.25

  def gen_dt(self, t):
    return t - self.t_n

  def gen_edot(self, e, t):
    return (e - self.e_n) / self.gen_dt(t)

  def gen_Tdot(self, T, t):
    return (T - self.T_n) / self.gen_dt(t)
