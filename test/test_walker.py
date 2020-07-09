import sys
sys.path.append('..')

from neml import walker, history, elasticity
from neml.math import tensors

from common import *
from nicediff import *
from test_visco_flow import CommonFlowRule
from test_general_flow import CommonGeneralFlow

import numpy.linalg as la

import unittest

class TestRateSwitch(unittest.TestCase, CommonGeneralFlow):
  def setUp(self):
    self.eps0 = 1.0e2
    self.D = 100.0
    self.n = 5.2
    self.s0 = 150.0
    self.K = 10000.0

    self.vmodel = walker.TestFlowRule(self.eps0, self.D, self.n, self.s0, self.K)

    E = 92000.0
    nu = 0.3

    self.emodel = elasticity.IsotropicLinearElasticModel(E, "youngs",
        nu, "poissons")
    
    self.lv = 1

    self.model = walker.WalkerKremplSwitchRule(self.emodel, self.vmodel,
        self.lv, self.eps0)

    self.T_n = 300.0
    self.e_n = np.zeros((6,))
    self.t_n = 0.0
    self.h_n = np.zeros((2,))

    self.eps = 1.0e-6

  def test_kappa(self):
    exact = self.model.kappa(self.gen_edot(self.gen_e(), self.gen_t()), self.gen_T())
    edev = make_dev(self.gen_edot(self.gen_e(), self.gen_t()))
    en = la.norm(edev)
    should = 1.0 - self.lv + self.lv * np.sqrt(2.0/3.0) * en / self.eps0
    self.assertAlmostEqual(should, exact)

  def test_dkappa(self):
    exact = self.model.dkappa(self.gen_edot(self.gen_e(), self.gen_t()), self.gen_T())
    should = differentiate(lambda e: self.model.kappa(e, self.gen_T()), 
        self.gen_edot(self.gen_e(), self.gen_t()))

    self.assertTrue(np.allclose(exact.flatten(),should.flatten(), rtol = 1e-4))

  def gen_hist(self):
    h = np.array([0.1,175.0])
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

class CommonWrappedFlow(object):
  def test_dy_ds(self):
    numerical = diff_scalar_symmetric(lambda s: self.model.y_wrap(self.make_state(s, self.h, self.T)), self.stress)
    actual = self.model.dy_ds_wrap(self.state)
    self.assertEqual(numerical, actual)

  def test_dy_dh(self):
    numerical = diff_history_scalar(lambda s: self.model.y_wrap(self.make_state(self.stress, s, self.T)), self.h)
    actual = self.model.dy_da_wrap(self.state)
    self.assertTrue(np.allclose(np.array(numerical), np.array(actual)))

  def test_dg_ds(self):
    numerical = diff_symmetric_symmetric(lambda s: self.model.g_wrap(self.make_state(s, self.h, self.T)), self.stress)
    actual = self.model.dg_ds_wrap(self.state)
    self.assertEqual(numerical, actual)

  def test_dg_dh(self):
    numerical = diff_symmetric_history(lambda h: self.model.g_wrap(self.make_state(self.stress, h, self.T)), self.h)
    actual= np.array(self.model.dg_da_wrap(self.state))
    self.assertTrue(np.allclose(np.array(numerical).reshape(actual.shape), actual))

  def test_dh_ds(self):
    numerical = diff_history_symmetric(lambda s: self.model.h_wrap(self.make_state(s, self.h, self.T)), self.stress)
    actual = np.array(self.model.dh_ds_wrap(self.state))
    self.assertTrue(np.allclose(np.array(numerical).reshape(actual.shape), actual))

  def test_dh_dh(self):
    numerical = diff_history_history(lambda h: self.model.h_wrap(self.make_state(self.stress, h, self.T)), self.h)
    actual = np.array(self.model.dh_da_wrap(self.state))
    self.assertTrue(np.allclose(np.array(numerical).reshape(actual.shape), actual))

  def test_dh_ds_time(self):
    numerical = diff_history_symmetric(lambda s: self.model.h_time_wrap(self.make_state(s, self.h, self.T)), self.stress)
    actual = np.array(self.model.dh_ds_time_wrap(self.state))
    self.assertTrue(np.allclose(np.array(numerical).reshape(actual.shape), actual))

  def test_dh_dh_time(self):
    numerical = diff_history_history(lambda h: self.model.h_time_wrap(self.make_state(self.stress, h, self.T)), self.h)
    actual = np.array(self.model.dh_da_time_wrap(self.state))
    self.assertTrue(np.allclose(np.array(numerical).reshape(actual.shape), actual))

  def test_dh_ds_temp(self):
    numerical = diff_history_symmetric(lambda s: self.model.h_temp_wrap(self.make_state(s, self.h, self.T)), self.stress)
    actual = np.array(self.model.dh_ds_temp_wrap(self.state))
    self.assertTrue(np.allclose(np.array(numerical).reshape(actual.shape), actual))

  def test_dh_dh_temp(self):
    numerical = diff_history_history(lambda h: self.model.h_temp_wrap(self.make_state(self.stress, h, self.T)), self.h)
    actual = np.array(self.model.dh_da_temp_wrap(self.state))
    self.assertTrue(np.allclose(np.array(numerical).reshape(actual.shape), actual))

class TestTestFlowRule(unittest.TestCase, CommonWrappedFlow, CommonFlowRule):
  def setUp(self):
    self.eps0 = 1.0e2
    self.D = 100.0
    self.n = 5.2
    self.s0 = 150.0
    self.K = 10000.0

    self.model = walker.TestFlowRule(self.eps0, self.D, self.n, self.s0, self.K)

    self.stress = tensors.Symmetric([
      [300.0,50.0,25.0],
      [50.0,150.0,-20.0],
      [25.0,-20.0,-100.0]])
    self.h = self.model.populate_hist()
    self.h.set_scalar("alpha", 0.1)
    self.h.set_scalar("iso", 200.0)
    self.T = 300.0

    self.state = self.make_state(self.stress, self.h, self.T)

    self.hist0 = np.array([0.0, self.s0])

  def gen_hist(self):
    return np.array([0.01,175.0])

  def make_state(self, S, h, T):
    return walker.State(S, h, T)

  def test_nhist(self):
    self.assertEqual(self.model.nhist, 2)

  def test_setup_hist(self):
    hobj = self.model.populate_hist()
    self.assertEqual(hobj.size, 2)
    self.assertTrue(hobj.items, ["alpha", "iso"])

  def test_initialize_hist(self):
    h = self.model.initialize_hist()
    self.assertAlmostEqual(h.get_scalar("alpha"), 0.0)
    self.assertAlmostEqual(h.get_scalar("iso"), self.s0)

  def test_y(self):
    should = self.eps0 * ((np.sqrt(3.0/2.0) * self.stress.dev().norm() - self.h.get_scalar("iso")) / self.D)**self.n
    actual = self.model.y_wrap(self.state)

    self.assertAlmostEqual(should, actual)

  def test_g(self):
    should = 3.0/2.0 * self.stress.dev() / (np.sqrt(3.0/2) * self.stress.dev().norm())
    actual = self.model.g_wrap(self.state)

    self.assertEqual(should, actual)

  def test_h(self):
    should = self.model.populate_hist()

    should.set_scalar("alpha", 1.0)
    should.set_scalar("iso", self.K)

    actual = self.model.h_wrap(self.state)

    self.assertTrue(np.allclose(np.array(should), np.array(actual)))

  def test_h_time(self):
    self.assertTrue(np.allclose(np.array(self.model.h_time_wrap(self.state)), np.zeros((2,))))

  def test_h_temp(self):
    self.assertTrue(np.allclose(np.array(self.model.h_temp_wrap(self.state)), np.zeros((2,))))

