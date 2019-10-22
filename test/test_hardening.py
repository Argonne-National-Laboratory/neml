import sys
sys.path.append('..')

from neml import hardening, interpolate
import unittest

from common import *

import numpy as np
import numpy.linalg as la

class CommonHardening(object):
  """
    Tests that can apply to all hardening rules
  """
  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def test_gradient(self):
    dfn = lambda x: self.model.q(x, self.T)
    ngrad = differentiate(dfn, self.hist_trial)
    grad = self.model.dq_da(self.hist_trial, self.T)
    self.assertTrue(np.allclose(ngrad, grad))

class TestLinearIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0

    self.hist0 = np.array([0.0])
    
    self.hist_trial = np.array([0.1])
    self.T = 300.0

    self.model = hardening.LinearIsotropicHardeningRule(self.s0, self.K)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.s0(self.T), self.s0))
    self.assertTrue(np.isclose(self.model.K(self.T), self.K))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      np.array([-self.s0 - self.K * self.hist_trial[0]])))

class TestInterpolatedIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.points = [0.0, 0.01, 0.02, 0.05]
    self.values = [100.0, 110.0, 130.0, 135.0]
    self.ifn = interpolate.PiecewiseLinearInterpolate(self.points, self.values)
    self.hist0 = np.array([0.0])
    self.hist_trial = np.array([0.03])
    self.T = 300.0

    self.model = hardening.InterpolatedIsotropicHardeningRule(self.ifn)

  def test_relation(self):
    self.assertTrue(np.isclose(self.model.q(self.hist_trial, self.T)[0],
      -self.ifn(self.hist_trial[0])))


class TestVoceIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.R = 100.0
    self.d = 10.0

    self.hist0 = np.array([0.0])
    
    self.hist_trial = np.abs([0.25])
    self.T = 300.0

    self.model = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.s0(self.T), self.s0))
    self.assertTrue(np.isclose(self.model.R(self.T), self.R))
    self.assertTrue(np.isclose(self.model.d(self.T), self.d))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      np.array([-self.s0 - self.R * (1 - np.exp(-self.d*self.hist_trial[0]))])))

class TestPowerLawIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.A = 100
    self.n = 0.2

    self.hist0 = np.array([0.0])
    
    self.hist_trial = np.abs([0.25])
    self.T = 300.0

    self.model = hardening.PowerLawIsotropicHardeningRule(self.s0, self.A, self.n)

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      np.array([-self.s0 - self.A * self.hist_trial[0]**self.n])))

class TestCombinedIsotropicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0

    self.hist0 = np.array([0.0])
    
    self.hist_trial = np.array([0.1])
    self.T = 300.0

    self.model1 = hardening.LinearIsotropicHardeningRule(self.s0, self.K)

    self.s0 = 200.0
    self.R = 100.0
    self.d = 10.0

    self.model2 = hardening.VoceIsotropicHardeningRule(self.s0, self.R, self.d)

    self.model = hardening.CombinedIsotropicHardeningRule([self.model1, self.model2])

  def test_combined(self):
    h1 = self.model1.q(self.hist_trial, self.T)
    h2 = self.model2.q(self.hist_trial, self.T)

    h = self.model.q(self.hist_trial, self.T)

    self.assertTrue(np.isclose(h, h1+h2))

class TestLinearKinematicHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.H = 1000.0

    self.hist0 = np.zeros((6,))
    
    self.hist_trial = np.array([50.0,25.0,75.0,10.0,15.0,50.0])
    self.hist_trial = self.hist_trial - np.array([1,1,1,0,0,0]) * sum(self.hist_trial[:3]) / 3.0
    self.T = 300.0

    self.model = hardening.LinearKinematicHardeningRule(self.H)

  def test_properties(self):
    self.assertTrue(np.isclose(self.model.H(self.T), self.H))

  def test_relation(self):
    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), 
      -self.hist_trial * self.H))

class TestCombinedHardening(unittest.TestCase, CommonHardening):
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0
    self.H = 1000.0

    self.hist0 = np.zeros((7,))
    self.hist_trial = np.array([50.0,75.0,100.0,25.0,10.0,15.0,60.0])
    self.hist_trial[1:] = make_dev(self.hist_trial[1:])
    self.T = 300.0

    self.iso = hardening.LinearIsotropicHardeningRule(self.s0, self.K)
    self.kin = hardening.LinearKinematicHardeningRule(self.H)

    self.model = hardening.CombinedHardeningRule(self.iso, self.kin)

  def test_relation(self):
    sb = np.zeros((7,))
    sb[0] = self.iso.q(self.hist_trial[0:1], self.T)
    sb[1:] = self.kin.q(self.hist_trial[1:], self.T)

    self.assertTrue(np.allclose(self.model.q(self.hist_trial, self.T), sb))

class CommonNonAssociative(object):
  """
    Common tests for non-associative hardening rules
  """
  def test_history(self):
    self.assertEqual(self.model.nhist, len(self.hist0))
    self.assertEqual(self.model.ninter, self.conform)
    self.assertTrue(np.allclose(self.model.init_hist(), self.hist0))

  def gen_stress(self):
    s = np.array([-150,200,-50,80,-20,30])
    return s

  def test_dq(self):
    a = self.gen_hist()

    dq_model = self.model.dq_da(a, self.T)
    dfn = lambda x: self.model.q(x, self.T)
    dq_num = differentiate(dfn, a)

    self.assertTrue(np.allclose(dq_model, dq_num))

  def test_dh_ds(self):
    s = self.gen_stress()
    a = self.gen_hist()

    dh_model = self.model.dh_ds(s, a, self.T)
    dfn = lambda x: self.model.h(x, a, self.T)
    dh_num = differentiate(dfn, s)

    self.assertTrue(np.allclose(dh_model, dh_num, rtol = 1.0e-3, atol = 1.0e-5))

  def test_dh_da(self):
    s = self.gen_stress()
    a = self.gen_hist()

    dh_model = self.model.dh_da(s, a, self.T)
    dfn = lambda x: self.model.h(s, x, self.T)
    dh_num = differentiate(dfn, a)

    self.assertTrue(np.allclose(dh_model, dh_num, rtol = 1.0e-3, atol = 1.0e-5))

  def test_dh_ds_time(self):
    s = self.gen_stress()
    a = self.gen_hist()

    dh_model = self.model.dh_ds_time(s, a, self.T)
    dfn = lambda x: self.model.h_time(x, a, self.T)
    dh_num = differentiate(dfn, s)

    self.assertTrue(np.allclose(dh_model, dh_num, rtol = 1.0e-3))

  def test_dh_da_time(self):
    s = self.gen_stress()
    a = self.gen_hist()

    dh_model = self.model.dh_da_time(s, a, self.T)
    dfn = lambda x: self.model.h_time(s, x, self.T)
    dh_num = differentiate(dfn, a)

    self.assertTrue(np.allclose(dh_model, dh_num, rtol = 1.0e-3, atol = 1.0e-4))

  def test_dh_ds_temp(self):
    s = self.gen_stress()
    a = self.gen_hist()

    dh_model = self.model.dh_ds_temp(s, a, self.T)
    dfn = lambda x: self.model.h_temp(x, a, self.T)
    dh_num = differentiate(dfn, s)

    self.assertTrue(np.allclose(dh_model, dh_num, rtol = 1.0e-3))

  def test_dh_da_temp(self):
    s = self.gen_stress()
    a = self.gen_hist()

    dh_model = self.model.dh_da_temp(s, a, self.T)
    dfn = lambda x: self.model.h_temp(s, x, self.T)
    dh_num = differentiate(dfn, a)

    self.assertTrue(np.allclose(dh_model, dh_num, rtol = 1.0e-3))

class CommonGamma(object):
  def gen_alpha(self):
    return 0.15

  def test_dgamma(self):
    a = self.gen_alpha()
    should = self.model.dgamma(a, self.T)
    dfn = lambda a: self.model.gamma(a, self.T)
    num = differentiate(dfn, a)

    self.assertTrue(np.isclose(should, num))

class TestConstantGamma(unittest.TestCase, CommonGamma):
  def setUp(self):
    self.g = 100.0
    self.model = hardening.ConstantGamma(self.g)

    self.T = 300.0

  def test_properties(self):
    self.assertTrue(np.isclose(self.g, self.model.g(self.T)))

  def test_gamma(self):
    a = self.gen_alpha()
    self.assertTrue(np.isclose(self.g, self.model.gamma(a, self.T)))

class TestSatGamma(unittest.TestCase, CommonGamma):
  def setUp(self):
    self.g0 = 100.0
    self.gs = 200.0
    self.beta = 2.5

    self.T = 300.0

    self.model = hardening.SatGamma(self.gs, self.g0, self.beta)

  def test_properties(self):
    self.assertTrue(np.isclose(self.g0, self.model.g0(self.T)))
    self.assertTrue(np.isclose(self.gs, self.model.gs(self.T)))
    self.assertTrue(np.isclose(self.beta, self.model.beta(self.T)))

  def test_gamma(self):
    a = self.gen_alpha()
    should = self.gs + (self.g0 - self.gs)*np.exp(-self.beta * a)
    self.assertTrue(np.isclose(should, self.model.gamma(a, self.T)))

class TestChaboche(unittest.TestCase, CommonNonAssociative):
  """
    Chaboche model with arbitrary hardening
  """
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0

    self.n = 4
    self.cs = np.array(range(self.n)) * 10.0
    self.rs = np.array(range(1,self.n+1)) * 10.0
    self.gammas = [hardening.ConstantGamma(r) for r in self.rs]

    self.As = (np.array(range(self.n)) + 1.0)  * 10.0
    self.a_s = np.array([2.0]*self.n)

    self.iso = hardening.LinearIsotropicHardeningRule(self.s0, self.K)

    self.model = hardening.Chaboche(self.iso, list(self.cs), self.gammas,
        list(self.As), list(self.a_s))

    self.hist0 = np.zeros((1 + self.n*6,))
    self.conform = 7

    self.T = 300.0

  def gen_hist(self):
    hist =  np.array(range(1,self.n*6+2)) / (self.n*7)
    hist[1:] = (1.0 - 2.0 * hist[1:]) * 100.0
    for i in range(self.n):
      hist[1+i*6:1+(i+1)*6] = make_dev(hist[1+i*6:1+(i+1)*6])
    return hist

  def test_properties(self):
    self.assertEqual(self.n, self.model.n)
    self.assertTrue(np.allclose(self.model.c(self.T), self.cs))

  def test_q(self):
    h = self.gen_hist()
    q_model = self.model.q(h, self.T)
    q_exact = np.zeros((7,))
    q_exact[0] = self.iso.q(h[0:1], self.T)
    q_exact[1:] = sum(h[1+i*6:1+(i+1)*6] for i in range(self.n))
    self.assertTrue(np.allclose(q_model, q_exact))

  def test_h(self):
    alpha = self.gen_hist()
    s = self.gen_stress()
    sdev = make_dev(s)
    X = sum(alpha[1+i*6:1+(i+1)*6] for i in range(self.n))
    n = (sdev+X) / la.norm(sdev+X)

    h_model = self.model.h(s, alpha, self.T)

    h_exact = np.zeros((self.model.nhist,))
    h_exact[0] = np.sqrt(2.0/3.0)
    for i in range(self.n):
      h_exact[1+i*6:1+(i+1)*6] = -2.0 / 3.0 * self.cs[i] * n - np.sqrt(2.0/3.0)*self.gammas[i].gamma(alpha[0], self.T)*alpha[1+i*6:1+(i+1)*6]
    
    self.assertTrue(np.allclose(h_model, h_exact))

  def test_h_time(self):
    alpha = self.gen_hist()
    s = self.gen_stress()
    sdev = make_dev(s)
    X = sum(alpha[1+i*6:1+(i+1)*6] for i in range(self.n))
    n = (sdev+X) / la.norm(sdev+X)

    h_model = self.model.h_time(s, alpha, self.T)

    h_exact = np.zeros((self.model.nhist,))
    for i in range(self.n):
      Xi = alpha[1+i*6:1+(i+1)*6] 
      h_exact[1+i*6:1+(i+1)*6] = -self.As[i] * np.sqrt(3.0/2.0) * la.norm(Xi)**(self.a_s[i]-1.0) * Xi
    
    self.assertTrue(np.allclose(h_model, h_exact))

class TestChabocheNewFormLinear(unittest.TestCase, CommonNonAssociative):
  """
    Same as previous but with new-style setup
  """
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0

    self.n = 4
    self.cs = list(np.array(range(self.n)) / (1+self.n))
    self.rs = np.array(range(self.n)) * 10.0
    self.gammas = [hardening.ConstantGamma(r) for r in self.rs]
    self.As = [0.0] * self.n
    self.ns = [1.0] * self.n

    self.iso = hardening.LinearIsotropicHardeningRule(self.s0, self.K)

    self.model = hardening.Chaboche(self.iso, self.cs, self.gammas, 
        self.As, self.ns)

    self.hist0 = np.zeros((1 + self.n*6,))
    self.conform = 7

    self.T = 300.0

  def gen_hist(self):
    hist = np.array(range(1,2+self.n*6)) / (3*self.n)
    hist[1:] = (1.0 - 2.0 * hist[1:]) * 100.0
    for i in range(self.n):
      hist[1+i*6:1+(i+1)*6] = make_dev(hist[1+i*6:1+(i+1)*6])
    return hist

  def test_properties(self):
    self.assertEqual(self.n, self.model.n)
    self.assertTrue(np.allclose(self.model.c(self.T), self.cs))

  def test_q(self):
    h = self.gen_hist()
    q_model = self.model.q(h, self.T)
    q_exact = np.zeros((7,))
    q_exact[0] = self.iso.q(h[0:1], self.T)
    q_exact[1:] = sum(h[1+i*6:1+(i+1)*6] for i in range(self.n))
    self.assertTrue(np.allclose(q_model, q_exact))

  def test_h(self):
    alpha = self.gen_hist()
    s = self.gen_stress()
    sdev = make_dev(s)
    X = sum(alpha[1+i*6:1+(i+1)*6] for i in range(self.n))
    n = (sdev+X) / la.norm(sdev+X)

    h_model = self.model.h(s, alpha, self.T)

    h_exact = np.zeros((self.model.nhist,))
    h_exact[0] = np.sqrt(2.0/3.0)
    for i in range(self.n):
      h_exact[1+i*6:1+(i+1)*6] = -2.0 / 3.0 * self.cs[i] * n - np.sqrt(2.0/3.0)*self.gammas[i].gamma(alpha[0], self.T)*alpha[1+i*6:1+(i+1)*6]
    
    self.assertTrue(np.allclose(h_model, h_exact))

class TestChabocheNewFormSat(unittest.TestCase, CommonNonAssociative):
  """
    Nonlinear saturation function
  """
  def setUp(self):
    self.s0 = 200.0
    self.K = 1000.0

    self.n = 2
    self.cs = list(np.array(range(self.n)) / (2*self.n)+10.0)
    self.g0s = np.array(range(1,self.n+1)) * 10.0
    self.gss = 2.0 * self.g0s
    self.betas = np.array(range(self.n)) / self.n
    self.gammas = [hardening.SatGamma(a, b, c) for a,b,c in zip(self.g0s, self.gss, self.betas)]
    self.As = [0.0] * self.n
    self.ns = [1.0] * self.n

    self.iso = hardening.LinearIsotropicHardeningRule(self.s0, self.K)

    self.model = hardening.Chaboche(self.iso, self.cs, self.gammas, self.As,
        self.ns)

    self.hist0 = np.zeros((1 + self.n*6,))
    self.conform = 7

    self.T = 300.0

  def gen_hist(self):
    hist = np.array(range(1, 2 + self.n*6))
    hist[1:] = (1.0 - 2.0 * hist[1:]) * 100.0
    for i in range(self.n):
      hist[1+i*6:1+(i+1)*6] = make_dev(hist[1+i*6:1+(i+1)*6])
    return hist

  def test_properties(self):
    self.assertEqual(self.n, self.model.n)
    self.assertTrue(np.allclose(self.model.c(self.T), self.cs))

  def test_q(self):
    h = self.gen_hist()
    q_model = self.model.q(h, self.T)
    q_exact = np.zeros((7,))
    q_exact[0] = self.iso.q(h[0:1], self.T)
    q_exact[1:] = sum(h[1+i*6:1+(i+1)*6] for i in range(self.n))
    self.assertTrue(np.allclose(q_model, q_exact))

  def test_h(self):
    alpha = self.gen_hist()
    s = self.gen_stress()
    sdev = make_dev(s)
    X = sum(alpha[1+i*6:1+(i+1)*6] for i in range(self.n))
    n = (sdev+X) / la.norm(sdev+X)

    h_model = self.model.h(s, alpha, self.T)

    h_exact = np.zeros((self.model.nhist,))
    h_exact[0] = np.sqrt(2.0/3.0)
    for i in range(self.n):
      h_exact[1+i*6:1+(i+1)*6] = -2.0 / 3.0 * self.cs[i] * n - np.sqrt(2.0/3.0)*self.gammas[i].gamma(alpha[0], self.T)*alpha[1+i*6:1+(i+1)*6]
    
    self.assertTrue(np.allclose(h_model, h_exact))

class TestChabocheTempTerm(unittest.TestCase, CommonNonAssociative):
  """
    Same as previous but with new-style setup
  """
  def setUp(self):
    s0 = 200.0
    K = 1000.0

    self.n = 2
    cvals1 = [1000.0,10.0]
    cvals2 = [100.0,100.0]
    Ts = [0.0,1000.0]

    cs = [interpolate.PiecewiseLinearInterpolate(Ts, [ci,cj]) for ci,cj in zip(cvals1, cvals2)]
    rs = [1.0e-2,1.0]
    gammas = [hardening.ConstantGamma(r) for r in rs]
    self.As = [0.0] * self.n
    self.ns = [1.0] * self.n

    iso = hardening.LinearIsotropicHardeningRule(s0, K)

    self.model = hardening.Chaboche(iso, cs, gammas, self.As, self.ns)

    self.hist0 = np.zeros((1 + self.n*6,))
    self.conform = 7

    self.T = 300.0
  
  def gen_hist(self):
    hist = np.array(range(1,2+self.n*6)) / (3*self.n)
    hist[1:] = (1.0 - 2.0 * hist[1:]) * 100.0
    for i in range(self.n):
      hist[1+i*6:1+(i+1)*6] = make_dev(hist[1+i*6:1+(i+1)*6])
    return hist
