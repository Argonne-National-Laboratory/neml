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

class CommonSofteningModel(object):
  def test_dphi(self):
    numerical = differentiate(lambda a: self.model.phi(a, self.T), self.a)
    actual = self.model.dphi(self.a, self.T)

    self.assertAlmostEqual(numerical, actual, delta = 1.0e-5)

class TestNoSoftening(unittest.TestCase, CommonSofteningModel):
  def setUp(self):
    self.a = 0.1
    self.T = 300.0

    self.model = walker.SofteningModel()

  def test_phi(self):
    self.assertAlmostEqual(self.model.phi(self.a, self.T), 1)

class TestWalkerSoftening(unittest.TestCase, CommonSofteningModel):
  def setUp(self):
    self.a = 0.1
    self.T = 300.0

    self.phi0 = 0.1
    self.phi1 = 2.1

    self.model = walker.WalkerSofteningModel(self.phi0, self.phi1)

  def test_phi(self):
    self.assertAlmostEqual(self.model.phi(self.a, self.T), 
        1.0 + self.phi0 * self.a**self.phi1)

class CommonThermalScaling(object):
  pass

class TestNoScaling(unittest.TestCase, CommonThermalScaling):
  def setUp(self):
    self.model = walker.ThermalScaling()
    self.T = 300.0

  def test_value(self):
    self.assertAlmostEqual(self.model.value(self.T), 1.0)

class TestArrheniusThermalScaling(unittest.TestCase, CommonThermalScaling):
  def setUp(self):
    self.Q = 64000.0
    self.R = 8.314
    self.T_ref = 300.0

    self.model = walker.ArrheniusThermalScaling(self.Q, self.R, self.T_ref)

    self.T = 523.0

  def test_value(self):
    should = np.exp(-self.Q/(self.R*self.T)) / np.exp(-self.Q/(self.R*self.T_ref))
    actual = self.model.value(self.T)
    self.assertAlmostEqual(should, actual)

class CommonIsotropicHardening(object):
  def make_state(self, h, a, adot, D, s, g, T):
    return walker.ScalarInternalVariableState(
        h, a, adot, D, s, g, T)

  def test_d_ratep_d_h(self):
    av = self.model.d_ratep_d_h(self.state)
    nv = differentiate(lambda h: 
        self.model.ratep(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertAlmostEqual(av, nv, places = 5)

  def test_d_ratep_d_a(self):
    av = self.model.d_ratep_d_a(self.state)
    nv = differentiate(lambda a: 
        self.model.ratep(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertAlmostEqual(av, nv)

  def test_d_ratep_d_adot(self):
    av = self.model.d_ratep_d_adot(self.state)
    nv = differentiate(lambda ad: 
        self.model.ratep(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_ratep_d_D(self):
    av = self.model.d_ratep_d_D(self.state)
    nv = differentiate(lambda D: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_ratep_d_s(self):
    av = self.model.d_ratep_d_s(self.state)
    nv = diff_scalar_symmetric(lambda s: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_ratep_d_g(self):
    av = self.model.d_ratep_d_g(self.state)
    nv = diff_scalar_symmetric(lambda g: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

  def test_ratet(self):
    self.assertAlmostEqual(self.model.ratet(self.state), 0)

  def test_d_ratet_d_h(self):
    av = self.model.d_ratet_d_h(self.state)
    nv = differentiate(lambda h: 
        self.model.ratet(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertAlmostEqual(av, nv, places = 5)

  def test_d_ratet_d_a(self):
    av = self.model.d_ratet_d_a(self.state)
    nv = differentiate(lambda a: 
        self.model.ratet(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertAlmostEqual(av, nv, places = 5)

  def test_d_ratet_d_adot(self):
    av = self.model.d_ratet_d_adot(self.state)
    nv = differentiate(lambda ad: 
        self.model.ratet(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv, places = 5)

  def test_d_ratet_d_D(self):
    av = self.model.d_ratet_d_D(self.state)
    nv = differentiate(lambda D: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_ratet_d_s(self):
    av = self.model.d_ratet_d_s(self.state)
    nv = diff_scalar_symmetric(lambda s: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_ratet_d_g(self):
    av = self.model.d_ratet_d_g(self.state)
    nv = diff_scalar_symmetric(lambda g: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

  def test_rateT(self):
    self.assertAlmostEqual(self.model.rateT(self.state), 0)

  def test_d_rateT_d_h(self):
    av = self.model.d_rateT_d_h(self.state)
    nv = differentiate(lambda h: 
        self.model.rateT(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertAlmostEqual(av, nv)

  def test_d_rateT_d_a(self):
    av = self.model.d_rateT_d_a(self.state)
    nv = differentiate(lambda a: 
        self.model.rateT(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertAlmostEqual(av, nv)

  def test_d_rateT_d_adot(self):
    av = self.model.d_rateT_d_adot(self.state)
    nv = differentiate(lambda ad: 
        self.model.rateT(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_ratet_d_D(self):
    av = self.model.d_rateT_d_D(self.state)
    nv = differentiate(lambda D: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_rateT_d_s(self):
    av = self.model.d_rateT_d_s(self.state)
    nv = diff_scalar_symmetric(lambda s: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_rateT_d_g(self):
    av = self.model.d_rateT_d_g(self.state)
    nv = diff_scalar_symmetric(lambda g: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

class TestConstantIsotropicHardening(unittest.TestCase,CommonIsotropicHardening):
  def setUp(self):
    self.model = walker.ConstantIsotropicHardening()

    self.h = 0.0
    self.a = 0.1
    self.adot = 2.0
    self.s = tensors.Symmetric([
      [300.0,50.0,25.0],
      [50.0,150.0,-20.0],
      [25.0,-20.0,-100.0]])
    self.g = self.s / self.s.norm()
    self.T = 350.0

    self.D = 10.0

    self.state = self.make_state(self.h, self.a, self.adot, self.D,
        self.s, self.g, self.T)

  def test_initial_value(self):
    self.assertAlmostEqual(self.model.initial_value(), 0.0)

  def test_ratep(self):
    self.assertAlmostEqual(self.model.ratep(self.state), 0)

class TestWalkerIsotropicHardening(unittest.TestCase,CommonIsotropicHardening):
  def setUp(self):
    self.r0 = 0.1
    self.r1 = 0.2
    self.r2 = 1.1
    self.Rinf = 0.25
    self.R0 = 0.01

    self.model = walker.WalkerIsotropicHardening(
        self.r0, self.Rinf, self.R0, self.r1, self.r2)

    self.h = 0.15
    self.a = 0.1
    self.adot = 2.0
    self.s = tensors.Symmetric([
      [300.0,50.0,25.0],
      [50.0,150.0,-20.0],
      [25.0,-20.0,-100.0]])
    self.g = self.s / self.s.norm()
    self.T = 350.0

    self.D = 10.0

    self.state = self.make_state(self.h, self.a, self.adot, self.D,
        self.s, self.g, self.T)

  def test_initial_value(self):
    self.assertAlmostEqual(self.model.initial_value(), 0.0)

  def test_ratep(self):
    self.assertAlmostEqual(self.model.ratep(self.state), 
        self.r0 * (self.Rinf - self.h))

  def test_ratet(self):
    self.assertAlmostEqual(self.model.ratet(self.state),
        self.r1 * (self.R0 - self.h) * np.abs(self.R0 - self.h)**(self.r2-1.0))

class CommonDragStress(object):
  def make_state(self, h, a, adot, D, s, g, T):
    return walker.ScalarInternalVariableState(
        h, a, adot, D, s, g, T)

  def test_d_ratep_d_h(self):
    av = self.model.d_ratep_d_h(self.state)
    nv = differentiate(lambda h: 
        self.model.ratep(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertAlmostEqual(av, nv, places = 5)

  def test_d_ratep_d_a(self):
    av = self.model.d_ratep_d_a(self.state)
    nv = differentiate(lambda a: 
        self.model.ratep(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertAlmostEqual(av, nv)

  def test_d_ratep_d_adot(self):
    av = self.model.d_ratep_d_adot(self.state)
    nv = differentiate(lambda ad: 
        self.model.ratep(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_ratep_d_D(self):
    av = self.model.d_ratep_d_D(self.state)
    nv = differentiate(lambda D: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_ratep_d_s(self):
    av = self.model.d_ratep_d_s(self.state)
    nv = diff_scalar_symmetric(lambda s: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_ratep_d_g(self):
    av = self.model.d_ratep_d_g(self.state)
    nv = diff_scalar_symmetric(lambda g: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

  def test_ratet(self):
    self.assertAlmostEqual(self.model.ratet(self.state), 0)

  def test_d_ratet_d_h(self):
    av = self.model.d_ratet_d_h(self.state)
    nv = differentiate(lambda h: 
        self.model.ratet(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertAlmostEqual(av, nv, places = 5)

  def test_d_ratet_d_a(self):
    av = self.model.d_ratet_d_a(self.state)
    nv = differentiate(lambda a: 
        self.model.ratet(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertAlmostEqual(av, nv, places = 5)

  def test_d_ratet_d_adot(self):
    av = self.model.d_ratet_d_adot(self.state)
    nv = differentiate(lambda ad: 
        self.model.ratet(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv, places = 5)

  def test_d_ratet_d_D(self):
    av = self.model.d_ratet_d_D(self.state)
    nv = differentiate(lambda D: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_ratet_d_s(self):
    av = self.model.d_ratet_d_s(self.state)
    nv = diff_scalar_symmetric(lambda s: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_ratet_d_g(self):
    av = self.model.d_ratet_d_g(self.state)
    nv = diff_scalar_symmetric(lambda g: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

  def test_rateT(self):
    self.assertAlmostEqual(self.model.rateT(self.state), 0)

  def test_d_rateT_d_h(self):
    av = self.model.d_rateT_d_h(self.state)
    nv = differentiate(lambda h: 
        self.model.rateT(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertAlmostEqual(av, nv)

  def test_d_rateT_d_a(self):
    av = self.model.d_rateT_d_a(self.state)
    nv = differentiate(lambda a: 
        self.model.rateT(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertAlmostEqual(av, nv)

  def test_d_rateT_d_adot(self):
    av = self.model.d_rateT_d_adot(self.state)
    nv = differentiate(lambda ad: 
        self.model.rateT(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_rateT_d_D(self):
    av = self.model.d_rateT_d_D(self.state)
    nv = differentiate(lambda D: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.adot)
    self.assertAlmostEqual(av, nv)

  def test_d_rateT_d_s(self):
    av = self.model.d_rateT_d_s(self.state)
    nv = diff_scalar_symmetric(lambda s: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_rateT_d_g(self):
    av = self.model.d_rateT_d_g(self.state)
    nv = diff_scalar_symmetric(lambda g: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

class TestConstantDragStress(CommonDragStress, unittest.TestCase):
  def setUp(self):
    self.value = 100.0
    self.model = walker.ConstantDragStress(self.value)

    self.h = 100.0
    self.a = 0.1
    self.adot = 2.0
    self.s = tensors.Symmetric([
      [300.0,50.0,25.0],
      [50.0,150.0,-20.0],
      [25.0,-20.0,-100.0]])
    self.g = self.s / self.s.norm()
    self.T = 350.0

    self.D = self.h

    self.state = self.make_state(self.h, self.a, self.adot, self.D,
        self.s, self.g, self.T)

  def test_initial_value(self):
    self.assertAlmostEqual(self.model.initial_value(), self.value)

  def test_ratep(self):
    self.assertAlmostEqual(self.model.ratep(self.state), 0)

  def test_D_xi(self):
    self.assertAlmostEqual(self.model.D_xi(self.T), 1.0)

  def test_D_0(self):
    self.assertAlmostEqual(self.model.D_0(self.T), self.value)

class TestWalkerDragStress(CommonDragStress, unittest.TestCase):
  def setUp(self):
    self.d0 = 4.856e6
    self.d1 = 8.55e-7
    self.d2 = 2.866
    self.D_xi = 121.6
    self.D_0 = 139.9 

    self.phi0 = 1.0    
    self.phi1 = 0.35

    self.Q = 3.65e5
    self.R = 8.314
    self.T0 = 950+273.15

    self.softening = walker.WalkerSofteningModel(self.phi0, self.phi1)
    self.scaling = walker.ArrheniusThermalScaling(self.Q, self.R, self.T0)

    self.model = walker.WalkerDragStress(self.d0,
        self.d1, self.d2, self.D_xi, self.D_0, 
        self.softening, scaling = self.scaling)

    self.h = 150.0
    self.a = 0.1
    self.adot = 2.0
    self.s = tensors.Symmetric([
      [300.0,50.0,25.0],
      [50.0,150.0,-20.0],
      [25.0,-20.0,-100.0]])
    self.g = self.s / self.s.norm()
    self.T = 900 + 273.15

    self.D = self.h

    self.state = self.make_state(self.h, self.a, self.adot, self.D,
        self.s, self.g, self.T)

  def test_initial_value(self):
    self.assertAlmostEqual(self.model.initial_value(), self.D_0)

  def test_ratep(self):
    should = self.d0 * (1.0 - (self.h-self.D_0) / self.D_xi)
    self.assertAlmostEqual(self.model.ratep(self.state), should)

  def test_ratet(self):
    should = -self.scaling.value(self.T) * self.softening.phi(self.a, self.T
        ) * self.d1 * (self.h-self.D_0)**(self.d2)
    self.assertAlmostEqual(self.model.ratet(self.state), should)

  def test_D_xi(self):
    self.assertAlmostEqual(self.model.D_xi(self.T), self.D_xi)

  def test_D_0(self):
    self.assertAlmostEqual(self.model.D_0(self.T), self.D_0)

class CommonKinematicHardening(object):
  def make_state(self, h, a, adot, D, s, g, T):
    return walker.SymmetricInternalVariableState(
        h, a, adot, D, s, g, T)

  def test_d_ratep_d_h(self):
    av = self.model.d_ratep_d_h(self.state)
    nv = diff_symmetric_symmetric(lambda h: 
        self.model.ratep(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertEqual(av, nv)

  def test_d_ratep_d_a(self):
    av = self.model.d_ratep_d_a(self.state)
    nv = diff_symmetric_scalar(lambda a: 
        self.model.ratep(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertEqual(av, nv)

  def test_d_ratep_d_adot(self):
    av = self.model.d_ratep_d_adot(self.state)
    nv =  diff_symmetric_scalar(lambda ad: 
        self.model.ratep(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertEqual(av, nv)

  def test_d_ratep_d_D(self):
    av = self.model.d_ratep_d_D(self.state)
    nv =  diff_symmetric_scalar(lambda D: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.adot)
    self.assertEqual(av, nv)

  def test_d_ratep_d_s(self):
    av = self.model.d_ratep_d_s(self.state)
    nv = diff_symmetric_symmetric(lambda s: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_ratep_d_g(self):
    av = self.model.d_ratep_d_g(self.state)
    nv = diff_symmetric_symmetric(lambda g: 
        self.model.ratep(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

  def test_ratet(self):
    self.assertEqual(self.model.ratet(self.state), tensors.Symmetric([[0,0,0],[0,0,0],[0,0,0]]))

  def test_d_ratet_d_h(self):
    av = self.model.d_ratet_d_h(self.state)
    nv = diff_symmetric_symmetric(lambda h: 
        self.model.ratet(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertEqual(av, nv)

  def test_d_ratet_d_a(self):
    av = self.model.d_ratet_d_a(self.state)
    nv =  diff_symmetric_scalar(lambda a: 
        self.model.ratet(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertEqual(av, nv)

  def test_d_ratet_d_adot(self):
    av = self.model.d_ratet_d_adot(self.state)
    nv =  diff_symmetric_scalar(lambda ad: 
        self.model.ratet(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertEqual(av, nv)

  def test_d_ratet_d_D(self):
    av = self.model.d_ratet_d_D(self.state)
    nv =  diff_symmetric_scalar(lambda D: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.D)
    self.assertEqual(av, nv)

  def test_d_ratet_d_s(self):
    av = self.model.d_ratet_d_s(self.state)
    nv = diff_symmetric_symmetric(lambda s: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_ratet_d_g(self):
    av = self.model.d_ratet_d_g(self.state)
    nv = diff_symmetric_symmetric(lambda g: 
        self.model.ratet(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

  def test_rateT(self):
    self.assertEqual(self.model.rateT(self.state), tensors.Symmetric([[0,0,0],[0,0,0],[0,0,0]]))

  def test_d_rateT_d_h(self):
    av = self.model.d_rateT_d_h(self.state)
    nv = diff_symmetric_symmetric(lambda h: 
        self.model.rateT(
          self.make_state(h, self.a, self.adot, self.D, self.s, self.g, self.T)),
        self.h)
    self.assertAlmostEqual(av, nv)

  def test_d_rateT_d_a(self):
    av = self.model.d_rateT_d_a(self.state)
    nv =  diff_symmetric_scalar(lambda a: 
        self.model.rateT(
          self.make_state(self.h, a, self.adot, self.D, self.s, self.g, self.T)),
        self.a)
    self.assertEqual(av, nv)

  def test_d_rateT_d_adot(self):
    av = self.model.d_rateT_d_adot(self.state)
    nv =  diff_symmetric_scalar(lambda ad: 
        self.model.rateT(
          self.make_state(self.h, self.a, ad, self.D, self.s, self.g, self.T)),
        self.adot)
    self.assertEqual(av, nv)

  def test_d_rateT_d_D(self):
    av = self.model.d_rateT_d_D(self.state)
    nv =  diff_symmetric_scalar(lambda D: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, D, self.s, self.g, self.T)),
        self.adot)
    self.assertEqual(av, nv)

  def test_d_rateT_d_s(self):
    av = self.model.d_rateT_d_s(self.state)
    nv = diff_symmetric_symmetric(lambda s: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, self.D, s, self.g, self.T)),
        self.s)
    self.assertEqual(av, nv)

  def test_d_rateT_d_g(self):
    av = self.model.d_rateT_d_g(self.state)
    nv = diff_symmetric_symmetric(lambda g: 
        self.model.rateT(
          self.make_state(self.h, self.a, self.adot, self.D, self.s, g, self.T)),
        self.g)
    self.assertEqual(av, nv)

class TestFAKinematicHardening(unittest.TestCase, CommonKinematicHardening):
  def setUp(self):
    self.c = 1000.0
    self.gamma = 2.1

    self.model = walker.FAKinematicHardening(self.c, self.gamma)

    self.h = tensors.Symmetric([
      [-100.0,200.0,10.0],
      [200.0,150.0,80.0],
      [10.0,80.0,-10.0]]).dev()
    self.a = 0.1
    self.adot = 2.0
    self.s = tensors.Symmetric([
      [300.0,50.0,25.0],
      [50.0,150.0,-20.0],
      [25.0,-20.0,-100.0]])
    self.g = self.s / self.s.norm()
    self.T = 900 + 273.15

    self.D = 10.0

    self.state = self.make_state(self.h, self.a, self.adot, self.D,
        self.s, self.g, self.T)

  def test_initial_value(self):
    self.assertEqual(self.model.initial_value(), tensors.Symmetric([[0,0,0],[0,0,0],[0,0,0]]))

  def test_ratep(self):
    self.assertEqual(self.model.ratep(self.state), 
        2.0/3.0 * self.c * self.g - self.gamma * self.h)

class TestWalkerKinematicHardening(unittest.TestCase, CommonKinematicHardening):
  def setUp(self):
    self.c0 = 1.723e5
    self.c1 = 3.153e6
    self.c2 = 1.691

    self.l0 = 14.0
    self.l1 = 0.7
    self.l = 6.165e1

    self.b0 = 0.5    # 0.0?
    self.x0 = 4.658e2
    self.x1 = 1.691e0

    self.phi0 = 0.5    # 0
    self.phi1 = 0.35

    self.Q = 3.65e5
    self.R = 8.314
    self.T0 = 950+273.15

    self.softening = walker.WalkerSofteningModel(self.phi0, self.phi1)
    self.scaling = walker.ArrheniusThermalScaling(self.Q, self.R, self.T0)

    self.model = walker.WalkerKinematicHardening(self.c0, self.c1, self.c2,
        self.l0, self.l1, self.l, self.b0, self.x0, self.x1, self.softening,
        scaling = self.scaling)

    self.h = tensors.Symmetric([
      [-100.0,200.0,10.0],
      [200.0,150.0,80.0],
      [10.0,80.0,-10.0]]).dev()
    self.a = 0.1
    self.adot = 2.0
    self.s = tensors.Symmetric([
      [300.0,50.0,25.0],
      [50.0,150.0,-20.0],
      [25.0,-20.0,-100.0]])
    self.g = self.s / self.s.norm()
    self.T = 900 + 273.15

    self.D = 75.0

    self.state = self.make_state(self.h, self.a, self.adot, self.D,
        self.s, self.g, self.T)

  def test_initial_value(self):
    self.assertEqual(self.model.initial_value(), tensors.Symmetric([[0,0,0],[0,0,0],[0,0,0]]))

  def test_ratep(self):
    c = self.c0 + self.c1 * self.adot**(1.0/self.c2)
    L = self.l * (self.l1 + (1.0-self.l1) * np.exp(-self.l0 * self.a))
    
    s = self.s.dev()
    n = 3.0/2.0 * (s - self.h) / (np.sqrt(3.0/2) * (s - self.h).norm())
    b = (1.0-self.b0) * self.h + 2.0/3.0 * self.b0 * tensors.douter(n, n).dot(self.h)

    self.assertEqual(self.model.ratep(self.state), 
        2.0/3.0 * c * self.g - c / L * b)

  def test_ratet(self):
    X = self.scaling.value(self.T)
    P = self.softening.phi(self.a, self.T)

    self.assertEqual(self.model.ratet(self.state),
        -X * self.x0 * P * (np.sqrt(3.0/2) * self.h.norm() / self.D)**self.x1 * self.h / (np.sqrt(3.0/2) * self.h.norm()))

class CommonWrappedFlow(object):
  def test_dy_ds_wrap(self):
    numerical = diff_scalar_symmetric(lambda s: self.model.y_wrap(self.make_state(s, self.h, self.T)), self.stress)
    actual = self.model.dy_ds_wrap(self.state)
    self.assertEqual(numerical, actual)

  def test_dy_dh_wrap(self):
    numerical = diff_history_scalar(lambda s: self.model.y_wrap(self.make_state(self.stress, s, self.T)), self.h)
    actual = self.model.dy_da_wrap(self.state)
      
    self.assertTrue(np.allclose(np.array(numerical), np.array(actual)))

  def test_dg_ds_wrap(self):
    numerical = diff_symmetric_symmetric(lambda s: self.model.g_wrap(self.make_state(s, self.h, self.T)), self.stress)
    actual = self.model.dg_ds_wrap(self.state)
    self.assertEqual(numerical, actual)

  def test_dg_dh_wrap(self):
    shape = (6, self.model.nhist)

    numerical = diff_symmetric_history(
        lambda h: self.model.g_wrap(self.make_state(self.stress, h, self.T)), self.h).reshape(shape)
    actual= np.array(self.model.dg_da_wrap(self.state)).reshape(shape[::-1]).T

    self.assertTrue(np.allclose(numerical, actual))

  def test_dh_ds_wrap(self):
    shape = (self.model.nhist, 6)

    numerical = diff_history_symmetric(lambda s: self.model.h_wrap(self.make_state(s, self.h, self.T)), self.stress).reshape(shape)
    actual = np.array(self.model.dh_ds_wrap(self.state)).reshape(shape)

    self.assertTrue(np.allclose(numerical, actual, rtol = 1e-3)) # HMM

  def test_dh_dh_wrap(self):
    shape = (self.model.nhist, self.model.nhist)

    numerical = diff_history_history(lambda h: self.model.h_wrap(self.make_state(self.stress, h, self.T)), self.h).reshape(shape)
    actual = self.model.dh_da_wrap(self.state).unravel_hh(self.h)
    
    self.assertTrue(np.allclose(numerical, actual, rtol = 1e-3))

  def test_dh_ds_time_wrap(self):
    shape = (self.model.nhist, 6)

    numerical = diff_history_symmetric(lambda s: self.model.h_time_wrap(self.make_state(s, self.h, self.T)), self.stress).reshape(shape)
    actual = np.array(self.model.dh_ds_time_wrap(self.state)).reshape(shape)

    self.assertTrue(np.allclose(numerical, actual, rtol = 1e-3)) # HMM

  def test_dh_dh_time_wrap(self):
    shape = (self.model.nhist, self.model.nhist)

    numerical = diff_history_history(lambda h: self.model.h_time_wrap(self.make_state(self.stress, h, self.T)), self.h).reshape(shape)
    actual = self.model.dh_da_time_wrap(self.state).unravel_hh(self.h)
    
    self.assertTrue(np.allclose(numerical, actual, rtol = 1e-3))

  def test_dh_ds_temp_wrap(self):
    shape = (self.model.nhist, 6)

    numerical = diff_history_symmetric(lambda s: self.model.h_temp_wrap(self.make_state(s, self.h, self.T)), self.stress).reshape(shape)
    actual = np.array(self.model.dh_ds_temp_wrap(self.state)).reshape(shape)

    self.assertTrue(np.allclose(numerical, actual, rtol = 1e-3)) # HMM

  def test_dh_dh_temp_wrap(self):
    shape = (self.model.nhist, self.model.nhist)

    numerical = diff_history_history(lambda h: self.model.h_temp_wrap(self.make_state(self.stress, h, self.T)), self.h).reshape(shape)
    actual = self.model.dh_da_temp_wrap(self.state).unravel_hh(self.h)
    
    self.assertTrue(np.allclose(numerical, actual, rtol = 1e-3))

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


class TestWalkerFlowRule(unittest.TestCase, CommonWrappedFlow, CommonFlowRule):
  def setUp(self):
    # Thermal scaling model
    self.Q = 3.65e5
    self.R = 8.314
    self.T0 = 950+273.15
    self.scaling = walker.ArrheniusThermalScaling(self.Q, self.R, self.T0)

    # Main softening model
    self.phi0 = 0.75  
    self.phi1 = 0.45
    self.Phi = walker.WalkerSofteningModel(self.phi0, self.phi1)

    # Reference strain rate
    self.eps0 = 1.0e-4

    # Rate sensitivity
    self.n = 6.1

    # Initial threshold
    self.k = 50.0

    # Switching function rate
    self.m = 1.2

    # Isotropic hardening
    self.r0 = 0.1
    self.r1 = 0.2
    self.r2 = 1.1
    self.Rinf = 0.25
    self.R0 = 0.01

    self.iso = walker.WalkerIsotropicHardening(
        self.r0, self.Rinf, self.R0, self.r1, self.r2)

    # Drag stress
    self.d0 = 4.856e6
    self.d1 = 8.55e-7
    self.d2 = 2.866
    self.D_xi = 121.6
    self.D_0 = 139.9 

    self.phi0_d = 1.0    
    self.phi1_d = 0.35

    self.softening_d = walker.WalkerSofteningModel(self.phi0_d, self.phi1_d)

    self.drag = walker.WalkerDragStress(self.d0,
        self.d1, self.d2, self.D_xi, self.D_0, 
        self.softening_d)

    # First back stress
    self.c01 = 1.723e5
    self.c11 = 3.153e6
    self.c21 = 1.691

    self.l01 = 14.0
    self.l11 = 0.7
    self.l1 = 6.165e1

    self.b01 = 0.5    # 0.0?
    self.x01 = 4.658e2
    self.x11 = 1.691e0

    self.phi01 = 0.5    # 0
    self.phi11 = 0.35

    self.softening1 = walker.WalkerSofteningModel(self.phi01, self.phi11)

    self.back1 = walker.WalkerKinematicHardening(self.c01, self.c11, self.c21,
        self.l01, self.l11, self.l1, self.b01, self.x01, self.x11, self.softening1)

    # Second back stress
    self.c02 = 2.723e5
    self.c12 = 1.153e6
    self.c22 = 2.691

    self.l02 = 12.0
    self.l12 = 0.8
    self.l2 = 4.165e1

    self.b02 = 0.6    # 0.0?
    self.x02 = 3.658e2
    self.x12 = 2.691e0

    self.phi02 = 0.8    # 0
    self.phi12 = 0.2

    self.softening2 = walker.WalkerSofteningModel(self.phi02, self.phi12)

    self.back2 = walker.WalkerKinematicHardening(self.c02, self.c12, self.c22,
        self.l02, self.l12, self.l2, self.b02, self.x02, self.x12, self.softening2)

    self.model = walker.WalkerFlowRule(self.eps0, self.Phi, self.scaling,
        self.n, self.k, self.m, self.iso, self.drag,  [self.back1, self.back2])

    self.stress = tensors.Symmetric([
      [300.0,50.0,25.0],
      [50.0,150.0,-20.0],
      [25.0,-20.0,-100.0]])
    
    self.a = 0.1
    self.R = 10.0
    self.D = 150.0

    self.X1 = tensors.Symmetric([
      [32.0,-10.0,18.0],
      [-10.0,25.0,-8.0],
      [18.0,-8.0,-10.0]]).dev()
    self.X2 = tensors.Symmetric([
      [-15.0,30.0,12.0],
      [30.0,-50.0,18.0],
      [12.0, 18.0,50.0]]).dev()

    self.X = self.X1 + self.X2

    self.h = self.model.populate_hist()
    self.h.set_scalar("alpha", self.a)
    self.h.set_scalar("R", self.R)
    self.h.set_scalar("D", self.D)
    self.h.set_symmetric("X0", self.X1)
    self.h.set_symmetric("X1", self.X2)
    self.T = 900.0+273.15

    self.state = self.make_state(self.stress, self.h, self.T)

    self.hist0 = np.array([0.0, 0.0, self.D_0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    self.szero = tensors.Symmetric([
      [0,0,0],
      [0.0,0,0],
      [0,0,0]])

    self.nhist = 15

  def gen_hist(self):
    return np.array([self.a, self.R, self.D, 9.07407407, 5.18518519, -14.25925926, -6.28539361, 14.14213562,
      -7.85674201, -8.33333333, -37.5, 45.83333333,  21.21320344, 14.14213562,  35.35533906])

  def make_state(self, S, h, T):
    return walker.State(S, h, T)

  def test_nhist(self):
    self.assertEqual(self.model.nhist, self.nhist)

  def test_setup_hist(self):
    hobj = self.model.populate_hist()
    self.assertEqual(hobj.size, self.nhist)
    self.assertTrue(hobj.items, ["alpha", "R", "D", "X0", "X1"])

  def test_initialize_hist(self):
    h = self.model.initialize_hist()
    self.assertAlmostEqual(h.get_scalar("alpha"), 0.0)
    self.assertAlmostEqual(h.get_scalar("R"), 0.0)
    self.assertAlmostEqual(h.get_scalar("D"), self.D_0)
    self.assertTrue(h.get_symmetric("X0"), self.szero)
    self.assertTrue(h.get_symmetric("X1"), self.szero)

  def test_y(self):
    mval = self.model.y_wrap(self.state)

    xi = (self.D - self.D_0) / self.D_xi
    Y = (self.k + self.R) * xi**self.m

    d = self.stress.dev() - self.X
    h = (np.sqrt(3.0/2) * d.norm() - Y) / self.D

    F = h**self.n

    a = self.eps0 * self.Phi.phi(self.a, self.T) * self.scaling.value(self.T)

    sval = a * F

    self.assertAlmostEqual(sval, mval)

  def test_g(self):
    mval = self.model.g_wrap(self.state)
    d = self.stress.dev() - self.X
    d = 3.0/2 * d / (np.sqrt(3.0/2.0) * d.norm())
    
    self.assertEqual(mval,d)

  def test_h(self):
    mval = self.model.h_wrap(self.state)
    sval = self.model.populate_hist()
    
    # All this for alpha...
    xi = (self.D - self.D_0) / self.D_xi
    Y = (self.k + self.R) * xi**self.m

    d = self.stress.dev() - self.X
    h = (np.sqrt(3.0/2) * d.norm() - Y) / self.D

    F = h**self.n

    a = self.eps0 * self.Phi.phi(self.a, self.T) * self.scaling.value(self.T)

    alpha = a * F
    
    sval.set_scalar("alpha", 1.0)

    # Need g....
    g = self.stress.dev() - self.X
    g = 3.0/2 * d / (np.sqrt(3.0/2.0) * d.norm())

    # R
    sval.set_scalar("R", self.iso.ratep(walker.ScalarInternalVariableState(
      self.R, self.a, alpha, self.D, self.stress, g, self.T)))

    # D
    sval.set_scalar("D", self.drag.ratep(walker.ScalarInternalVariableState(
      self.D, self.a, alpha, self.D, self.stress, g, self.T)))

    # X1
    sval.set_symmetric("X0", self.back1.ratep(walker.SymmetricInternalVariableState(
      self.X1, self.a, alpha, self.D, self.stress, g, self.T)))

    # X2
    sval.set_symmetric("X1", self.back2.ratep(walker.SymmetricInternalVariableState(
      self.X2, self.a, alpha, self.D, self.stress, g, self.T)))
    
    self.assertTrue(np.allclose(np.array(sval), np.array(mval)))

  def test_h_time(self):
    mval = self.model.h_time_wrap(self.state)
    sval = self.model.populate_hist()
    
    # All this for alpha...
    xi = (self.D - self.D_0) / self.D_xi
    Y = (self.k + self.R) * xi**self.m

    d = self.stress.dev() - self.X
    h = (np.sqrt(3.0/2) * d.norm() - Y) / self.D

    F = h**self.n

    a = self.eps0 * self.Phi.phi(self.a, self.T) * self.scaling.value(self.T)

    alpha = a * F
    
    # Time rate...
    sval.set_scalar("alpha", 0.0)

    # Need g....
    g = self.stress.dev() - self.X
    g = 3.0/2 * d / (np.sqrt(3.0/2.0) * d.norm())

    # R
    sval.set_scalar("R", self.iso.ratet(walker.ScalarInternalVariableState(
      self.R, self.a, alpha, self.D, self.stress, g, self.T)))

    # D
    sval.set_scalar("D", self.drag.ratet(walker.ScalarInternalVariableState(
      self.D, self.a, alpha, self.D, self.stress, g, self.T)))

    # X1
    sval.set_symmetric("X0", self.back1.ratet(walker.SymmetricInternalVariableState(
      self.X1, self.a, alpha, self.D, self.stress, g, self.T)))

    # X2
    sval.set_symmetric("X1", self.back2.ratet(walker.SymmetricInternalVariableState(
      self.X2, self.a, alpha, self.D, self.stress, g, self.T)))
    
    self.assertTrue(np.allclose(np.array(sval), np.array(mval)))

  def test_h_temp(self):
    self.assertTrue(np.allclose(np.array(self.model.h_temp_wrap(self.state)), np.zeros((self.nhist,))))

