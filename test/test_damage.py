import sys
sys.path.append('..')

from neml import interpolate, solvers, models, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep, damage, larsonmiller
from common import *

import unittest
import numpy as np
import numpy.linalg as la

class CommonStandardDamageModel(object):
  """
    Tests that apply to any standard damage model
  """
  def effective(self, s):
    sdev = make_dev(s)
    return np.sqrt(3.0/2.0 * np.dot(sdev, sdev))

  def test_damage(self):
    d_model = self.model.damage(self.d_np1, self.d_n, self.e_np1, self.e_n,
        self.s_np1, self.s_n, self.T_np1, self.T_n, self.t_np1, self.t_n)
    S = self.elastic.S(self.T_np1)
    dS = self.s_np1 - self.s_n
    dee = np.dot(S, dS)
    de = self.e_np1 - self.e_n
    
    dp = np.sqrt(2.0/3.0 * (np.dot(de, de) + np.dot(dee, dee) - 
      2.0 * np.dot(dee, de)))
    f = self.model.f(self.s_np1, self.d_np1, self.T_np1)
    d_calcd = self.d_n + f * dp

    self.assertTrue(np.isclose(d_model, d_calcd))

  def test_function_derivative_s(self):
    d_model = self.model.df_ds(self.stress, self.d_np1, self.T)
    d_calcd = differentiate(lambda x: self.model.f(x, self.d_np1, self.T),
        self.stress)
    self.assertTrue(np.allclose(d_model, d_calcd))

  def test_function_derivative_d(self):
    d_model = self.model.df_dd(self.stress, self.d, self.T)
    d_calcd = differentiate(lambda x: self.model.f(self.s_np1, x, self.T),
        self.d)
    self.assertTrue(np.isclose(d_model, d_calcd))

class CommonScalarDamageModel(object):
  def test_ndamage(self):
    self.assertEqual(self.model.ndamage, 1)

  def test_init_damage(self):
    self.assertTrue(np.allclose(self.model.init_damage(), np.zeros((1,))))

  def test_ddamage_ddamage(self):
    dd_model = self.model.ddamage_dd(self.d_np1, self.d_n, self.e_np1, self.e_n,
        self.s_np1, self.s_n, self.T_np1, self.T_n, self.t_np1, self.t_n)
    dfn = lambda d: self.model.damage(d, self.d_n, self.e_np1, self.e_n,
        self.s_np1, self.s_n, self.T_np1, self.T_n, self.t_np1, self.t_n)
    dd_calcd = differentiate(dfn, self.d_np1)

    print(dd_model, dd_calcd)
    
    self.assertTrue(np.isclose(dd_model, dd_calcd))

  def test_ddamage_dstrain(self):
    dd_model = self.model.ddamage_de(self.d_np1, self.d_n, self.e_np1, self.e_n,
        self.s_np1, self.s_n, self.T_np1, self.T_n, self.t_np1, self.t_n)
    dfn = lambda e: self.model.damage(self.d_np1, self.d_n, e, self.e_n,
        self.s_np1, self.s_n, self.T_np1, self.T_n, self.t_np1, self.t_n)
    dd_calcd = differentiate(dfn, self.e_np1)[0]

    self.assertTrue(np.allclose(dd_model, dd_calcd, rtol = 1.0e-3))

  def test_ddamage_dstress(self):
    dd_model = self.model.ddamage_ds(self.d_np1, self.d_n, self.e_np1, self.e_n,
        self.s_np1, self.s_n, self.T_np1, self.T_n, self.t_np1, self.t_n)
    dfn = lambda s: self.model.damage(self.d_np1, self.d_n, self.e_np1, self.e_n,
        s, self.s_n, self.T_np1, self.T_n, self.t_np1, self.t_n)
    dd_calcd = differentiate(dfn, self.s_np1)[0]

    self.assertTrue(np.allclose(dd_model, dd_calcd, rtol = 1.0e-3))

  def test_nparams(self):
    self.assertEqual(self.model.nparams, 7)

  def test_init_x(self):
    trial_state = self.model.make_trial_state(
        self.e_np1, self.e_n,
        self.T_np1, self.T_n, self.t_np1, self.t_n,
        self.s_n, self.hist_n, self.u_n, self.p_n)
    
    me = np.array(list(self.s_n) + [self.d_n])
    them = self.model.init_x(trial_state)

    self.assertTrue(np.allclose(me, them))

  def test_R(self):
    trial_state = self.model.make_trial_state(
        self.e_np1, self.e_n,
        self.T_np1, self.T_n, self.t_np1, self.t_n,
        self.s_n, self.hist_n, self.u_n, self.p_n)

    R, J = self.model.RJ(self.x_trial, trial_state)

    s_trial = self.x_trial[:6]
    w_trial = self.x_trial[6]

    R_calc = np.zeros((7,))
    s_p_np1, h_p, A_p, u_p, p_p =self.bmodel.update_sd(self.e_np1, self.e_n,
        self.T_np1, self.T_n,
        self.t_np1, self.t_n, self.s_n / (1-self.d_n), self.hist_n[1:],
        self.u_n, self.p_n)
    R_calc[:6] = s_trial - (1-w_trial) * s_p_np1
    d_np1 = self.model.damage(w_trial, self.d_n, self.e_np1, self.e_n, 
        s_trial / (1 - w_trial), self.s_n / (1-self.d_n), self.T_np1,
        self.T_n, self.t_np1, self.t_n)
    R_calc[6] = w_trial - d_np1

    self.assertTrue(np.allclose(R_calc, R))

  def test_jacobian(self):
    trial_state = self.model.make_trial_state(
        self.e_np1, self.e_n,
        self.T_np1, self.T_n, self.t_np1, self.t_n,
        self.s_n, self.hist_n, self.u_n, self.p_n)

    R, J = self.model.RJ(self.x_trial, trial_state)
    Jnum = differentiate(lambda x: self.model.RJ(x, trial_state)[0], 
        self.x_trial)
    
    self.assertTrue(np.allclose(J, Jnum, rtol = 1.0e-3))

class CommonDamagedModel(object):
  def test_nstore(self):
    self.assertEqual(self.model.nstore, self.bmodel.nstore + self.model.ndamage)

  def test_store(self):
    base = self.bmodel.init_store()
    damg = self.model.init_damage()
    comp = list(damg) + list(base)
    fromm = self.model.init_store()
    self.assertTrue(np.allclose(fromm, comp))
  
  def test_tangent_proportional_strain(self):
    t_n = 0.0
    e_n = np.zeros((6,))
    s_n = np.zeros((6,))
    hist_n = self.model.init_store()
    u_n = 0.0
    p_n = 0.0

    for m in np.linspace(0,1,self.nsteps+1)[1:]:
      t_np1 = m * self.ttarget
      e_np1 = m * self.etarget
     
      trial_state = self.model.make_trial_state(
          e_np1, e_n,
          self.T, self.T, t_np1, t_n,
          s_n, hist_n, u_n, p_n)

      s_np1, hist_np1, A_np1, u_np1, p_np1 = self.model.update_sd(
          e_np1, e_n, self.T, self.T, t_np1, t_n, s_n, hist_n,
          u_n, p_n)

      A_num = differentiate(lambda e: self.model.update_sd(e, e_n,
        self.T, self.T, t_np1, t_n, s_n, hist_n, u_n, p_n)[0], e_np1)
      
      print(A_num)
      print(A_np1)
      self.assertTrue(np.allclose(A_num, A_np1, rtol = 5.0e-2, atol = 1.0e-1))

      e_n = np.copy(e_np1)
      s_n = np.copy(s_np1)
      hist_n = np.copy(hist_np1)
      u_n = u_np1
      p_n = p_np1
      t_n = t_np1

class TestWorkDamage(unittest.TestCase, CommonScalarDamageModel,
    CommonDamagedModel):
  def setUp(self):
    self.E = 92000.0
    self.nu = 0.3

    self.s0 = 180.0
    self.Kp = 1000.0
    self.H = 1000.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E, "youngs",
        self.nu, "poissons")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.bmodel = models.SmallStrainRateIndependentPlasticity(self.elastic, 
        flow)
    
    self.fn = interpolate.PolynomialInterpolate([0.1,5.0, 1e-8])
    self.n = 2.1

    self.model = damage.NEMLWorkDamagedModel_sd(
        self.elastic, self.fn, self.n, self.bmodel, verbose = False)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0]) * 0.75
    self.T = 100.0

    self.s_np1 = self.stress
    self.s_n = np.array([-25,150,250,-25,-100,25]) * 0.0

    self.d_np1 = 0.2
    self.d_n = 0.1

    self.e_np1 = np.array([0.1,-0.01,0.15,-0.05,-0.1,0.15]) * 0.75
    self.e_n = np.array([-0.05,0.025,-0.1,0.2,0.11,0.13]) * 0.0

    self.T_np1 = self.T
    self.T_n = 90.0

    self.t_np1 = 5000.0
    self.t_n = 0.9

    self.u_n = 0.0
    self.p_n = 0.0
  
    # This is a rather boring baseline history state to probe, but I can't
    # think of a better way to get a "generic" history from a generic model
    self.hist_n = np.array([self.d_n] + list(self.bmodel.init_store()))
    self.x_trial = np.array([50,-25,150,-150,190,100.0] + [0.41])

    self.nsteps = 10
    self.etarget = np.array([0.1,-0.025,0.02,0.015,-0.02,-0.05])
    self.ttarget = 2.0

    self.ee = np.dot(self.elastic.S(self.T_np1), self.s_np1*(1.0-self.d_np1) - self.s_n*(1.0-self.d_n))
    self.de = self.e_np1 - self.e_n
    self.dp = self.de - self.ee
    self.dt = self.t_np1 - self.t_n
    self.Wdot = np.dot(self.s_np1*(1.0-self.d_np1), self.dp) / self.dt

  def test_definition(self):
    damage = self.model.damage(self.d_np1, self.d_n, self.e_np1, self.e_n,
        self.s_np1, self.s_n, self.T_np1, self.T_n, self.t_np1, self.t_n)
    should = self.d_n + self.n * self.d_np1**((self.n-1)/self.n) * self.Wdot * self.dt / self.fn.value(self.Wdot) 

    self.assertAlmostEqual(damage, should)

class TestClassicalDamage(unittest.TestCase, CommonScalarDamageModel,
    CommonDamagedModel):
  def setUp(self):
    self.E = 92000.0
    self.nu = 0.3

    self.s0 = 180.0
    self.Kp = 1000.0
    self.H = 1000.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E, "youngs",
        self.nu, "poissons")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.bmodel = models.SmallStrainRateIndependentPlasticity(self.elastic, 
        flow)

    self.xi = 0.478
    self.phi = 1.914
    self.A = 10000000.0

    self.model = damage.ClassicalCreepDamageModel_sd(
        self.elastic, 
        self.A, self.xi, self.phi, self.bmodel)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])
    self.T = 100.0

    self.s_np1 = self.stress
    self.s_n = np.array([-25,150,250,-25,-100,25])

    self.d_np1 = 0.5
    self.d_n = 0.4

    self.e_np1 = np.array([0.1,-0.01,0.15,-0.05,-0.1,0.15])
    self.e_n = np.array([-0.05,0.025,-0.1,0.2,0.11,0.13])

    self.T_np1 = self.T
    self.T_n = 90.0

    self.t_np1 = 1.0
    self.t_n = 0.0

    self.u_n = 0.0
    self.p_n = 0.0
  
    # This is a rather boring baseline history state to probe, but I can't
    # think of a better way to get a "generic" history from a generic model
    self.hist_n = np.array([self.d_n] + list(self.bmodel.init_store()))
    self.x_trial = np.array([50,-25,150,-150,190,100.0] + [0.41])

    self.nsteps = 10
    self.etarget = np.array([0.1,-0.025,0.02,0.015,-0.02,-0.05])
    self.ttarget = 10.0

class CommonEffectiveStress(object):
  def test_derivative(self):
    nd = differentiate(self.es.effective, self.stress)
    d = self.es.deffective(self.stress)
    self.assertTrue(np.allclose(nd, d))

class TestVonMisesEffectiveStress(unittest.TestCase, CommonEffectiveStress):
  def setUp(self):
    self.es = damage.VonMisesEffectiveStress()

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])

  def test_definition(self):
    S = usym(self.stress)
    D = S - np.trace(S) * np.eye(3) / 3.0
    se = np.sqrt(3.0/2 * np.sum(D*D))
    self.assertTrue(np.isclose(se, self.es.effective(self.stress)))

class TestHuddlestonEffectiveStress(unittest.TestCase, CommonEffectiveStress):
  def setUp(self):
    self.b = 0.25
    self.es = damage.HuddlestonEffectiveStress(self.b)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])

  def test_definition(self):
    S = usym(self.stress)
    s3, s2, s1 = la.eigvalsh(S)

    J1 = s1 + s2 + s3
    sb = np.sqrt(((s1-s2)**2.0+(s2-s3)**2.0+(s3-s1)**2.0)/2.0)
    ss = np.sqrt(s1**2.0 + s2**2.0 + s3**2.0)

    v1 = sb * np.exp(self.b*(J1/ss-1.0))

    I1 = np.trace(S)
    I2 = 0.5*(np.trace(S)**2.0 - np.trace(np.dot(S,S)))
    Sd = S - np.trace(S) * np.eye(3) / 3.0
    I2p = 0.5*(np.trace(Sd)**2.0 - np.trace(np.dot(Sd,Sd)))
    
    J2p = -I2p # Very odd thing here

    sb2 = np.sqrt(3.0*J2p)
    ss2 = np.sqrt(3.0*J2p+I2)

    v2 = sb2 * np.exp(self.b*(I1/ss2-1.0))
    v3 = self.es.effective(self.stress)
    
    self.assertTrue(np.isclose(sb,sb2))
    self.assertTrue(np.isclose(ss,ss2))
    
    self.assertTrue(np.isclose(v1,v2))
    self.assertTrue(np.isclose(v1,v3))

class TestMaxPrincipalEffectiveStress(unittest.TestCase, CommonEffectiveStress):
  def setUp(self):
    self.es = damage.MaxPrincipalEffectiveStress()

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])

  def test_definition(self):
    S = usym(self.stress)
    vals = la.eigvalsh(S)
    se = np.max(vals)
    self.assertTrue(np.isclose(se, self.es.effective(self.stress)))

class TestMaxSeveralEffectiveStress(unittest.TestCase, CommonEffectiveStress):
  def setUp(self):
    self.es = damage.MaxSeveralEffectiveStress([
      damage.VonMisesEffectiveStress(),
      damage.MaxPrincipalEffectiveStress()])

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])

  def test_definition(self):
    S = usym(self.stress)
    vals = la.eigvalsh(S)
    se1 = np.max(vals)

    D = S - np.trace(S) * np.eye(3) / 3.0
    se2 = np.sqrt(3.0/2 * np.sum(D*D))

    self.assertTrue(np.isclose(max(se1,se2), self.es.effective(self.stress)))

class TestSumSeveralEffectiveStress(unittest.TestCase, CommonEffectiveStress):
  def setUp(self):
    self.stresses = [damage.VonMisesEffectiveStress(), 
        damage.MaxPrincipalEffectiveStress()]
    self.weights = [0.37, 1.0-0.37]
    self.es = damage.SumSeveralEffectiveStress(self.stresses, self.weights)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])

  def test_definition(self):
    base = sum(w * eff.effective(self.stress) for eff, w in zip(self.stresses,
      self.weights))

    self.assertTrue(np.isclose(base, self.es.effective(self.stress)))

class BaseModularDamage(CommonScalarDamageModel, CommonDamagedModel):
  def complete(self):
    self.E = 92000.0
    self.nu = 0.3

    self.s0 = 180.0
    self.Kp = 1000.0
    self.H = 1000.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E, "youngs",
        self.nu, "poissons")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.bmodel = models.SmallStrainRateIndependentPlasticity(self.elastic, 
        flow)

    self.xi = 0.478
    self.phi = 1.914
    self.A = 10000000.0

    self.model = damage.ModularCreepDamageModel_sd(
        self.elastic, 
        self.A, self.xi, self.phi,
        self.effective_model(),
        self.bmodel)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])
    self.T = 100.0

    self.s_np1 = self.stress
    self.s_n = np.array([-25,150,250,-25,-100,25])

    self.d_np1 = 0.5
    self.d_n = 0.4

    self.e_np1 = np.array([0.1,-0.01,0.15,-0.05,-0.1,0.15])
    self.e_n = np.array([-0.05,0.025,-0.1,0.2,0.11,0.13])

    self.T_np1 = self.T
    self.T_n = 90.0

    self.t_np1 = 1.0
    self.t_n = 0.0

    self.u_n = 0.0
    self.p_n = 0.0
  
    # This is a rather boring baseline history state to probe, but I can't
    # think of a better way to get a "generic" history from a generic model
    self.hist_n = np.array([self.d_n] + list(self.bmodel.init_store()))
    self.x_trial = np.array([50,-25,150,-150,190,100.0] + [0.41])

    self.nsteps = 10
    self.etarget = np.array([0.1,-0.025,0.02,0.015,-0.02,-0.05])
    self.ttarget = 10.0

class TestModularVonMises(unittest.TestCase, BaseModularDamage):
  def setUp(self):
    self.complete()

  def effective_model(self):
    return damage.VonMisesEffectiveStress()

class TestModularMaxPrincipal(unittest.TestCase, BaseModularDamage):
  def setUp(self):
    self.complete()

  def effective_model(self):
    return damage.MaxPrincipalEffectiveStress()

class TestModularMaxSeveral(unittest.TestCase, BaseModularDamage):
  def setUp(self):
    self.complete()

  def effective_model(self):
    return damage.MaxSeveralEffectiveStress([
      damage.VonMisesEffectiveStress(),
      damage.MaxPrincipalEffectiveStress()])

class TestHuddlestonPrincipal(unittest.TestCase, BaseModularDamage):
  def setUp(self):
    self.complete()

  def effective_model(self):
    return damage.HuddlestonEffectiveStress(0.24)

class TestLMDamage(unittest.TestCase, CommonScalarDamageModel, CommonDamagedModel):
  def setUp(self):
    self.E = 92000.0
    self.nu = 0.3

    self.s0 = 180.0
    self.Kp = 1000.0
    self.H = 1000.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E, "youngs",
        self.nu, "poissons")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.bmodel = models.SmallStrainRateIndependentPlasticity(self.elastic, 
        flow)

    self.fn = interpolate.PolynomialInterpolate([-6.653e-9,2.952e-4,-6.197e-1])
    self.C = 32.06

    self.lmr = larsonmiller.LarsonMillerRelation(self.fn, self.C)
    self.effective = damage.VonMisesEffectiveStress()

    self.model = damage.LarsonMillerCreepDamageModel_sd(
        self.elastic, self.lmr, self.effective, self.bmodel)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])
    self.T = 100.0

    self.s_np1 = self.stress
    self.s_n = np.array([-25,150,250,-25,-100,25])

    self.d_np1 = 0.5
    self.d_n = 0.4

    self.e_np1 = np.array([0.1,-0.01,0.15,-0.05,-0.1,0.15])
    self.e_n = np.array([-0.05,0.025,-0.1,0.2,0.11,0.13])

    self.T_np1 = self.T
    self.T_n = 90.0

    self.t_np1 = 1.0
    self.t_n = 0.0
    self.dt = self.t_np1 - self.t_n

    self.u_n = 0.0
    self.p_n = 0.0
  
    # This is a rather boring baseline history state to probe, but I can't
    # think of a better way to get a "generic" history from a generic model
    self.hist_n = np.array([self.d_n] + list(self.bmodel.init_store()))
    self.x_trial = np.array([50,-25,150,-150,190,100.0] + [0.41])

    self.nsteps = 10
    self.etarget = np.array([0.1,-0.025,0.02,0.015,-0.02,-0.05])
    self.ttarget = 10.0

  def test_damage(self):
    se = self.effective.effective(self.stress)
    tR = self.lmr.tR(se, self.T)
    exact = self.d_n + self.dt / tR
    model = self.model.damage(self.d_np1, self.d_n, 
        self.e_np1, self.e_n, self.s_np1, self.s_n,
        self.T_np1, self.T_n, self.t_np1, self.t_n)
    self.assertTrue(np.isclose(exact, model))

class TestPowerLawDamage(unittest.TestCase, CommonStandardDamageModel, 
    CommonScalarDamageModel, CommonDamagedModel):
  def setUp(self):
    self.E = 92000.0
    self.nu = 0.3

    self.s0 = 180.0
    self.Kp = 1000.0
    self.H = 1000.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E, "youngs",
        self.nu, "poissons")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.bmodel = models.SmallStrainRateIndependentPlasticity(self.elastic, 
        flow)

    self.A = 8.0e-6
    self.a = 2.2

    self.model = damage.NEMLPowerLawDamagedModel_sd(self.elastic, self.A, self.a, 
        self.bmodel)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])
    self.T = 100.0
    self.d = 0.45

    self.s_np1 = self.stress
    self.s_n = np.array([-25,150,250,-25,-100,25])

    self.d_np1 = 0.5
    self.d_n = 0.4

    self.e_np1 = np.array([0.1,-0.01,0.15,-0.05,-0.1,0.15])
    self.e_n = np.array([-0.05,0.025,-0.1,0.2,0.11,0.13])

    self.T_np1 = self.T
    self.T_n = 90.0

    self.t_np1 = 1.0
    self.t_n = 0.0

    self.u_n = 0.0
    self.p_n = 0.0
  
    # This is a rather boring baseline history state to probe, but I can't
    # think of a better way to get a "generic" history from a generic model
    self.hist_n = np.array([self.d_n] + list(self.bmodel.init_store()))
    self.x_trial = np.array([50,-25,150,-150,190,100.0] + [0.41])

    self.nsteps = 10
    self.etarget = np.array([0.1,-0.025,0.02,0.015,-0.02,-0.05])
    self.ttarget = 10.0

  def test_function(self):
    f_model = self.model.f(self.stress, self.d_np1, self.T)
    f_calcd = self.A * self.effective(self.stress) ** self.a
    self.assertTrue(np.isclose(f_model, f_calcd))

class TestExponentialDamage(unittest.TestCase, CommonStandardDamageModel, 
    CommonScalarDamageModel, CommonDamagedModel):
  def setUp(self):
    self.E = 92000.0
    self.nu = 0.3

    self.s0 = 180.0
    self.Kp = 1000.0
    self.H = 1000.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E,
        "youngs", self.nu, "poissons")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.bmodel = models.SmallStrainRateIndependentPlasticity(self.elastic, 
        flow)

    self.W0 = 10.0
    self.k0 = 0.0001
    self.a = 2.0

    self.model = damage.NEMLExponentialWorkDamagedModel_sd(
        self.elastic, self.W0, self.k0, 
        self.a, self.bmodel)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])
    self.T = 100.0
    self.d = 0.45

    self.s_np1 = self.stress
    self.s_n = np.array([-25,150,250,-25,-100,25])

    self.d_np1 = 0.5
    self.d_n = 0.4

    self.e_np1 = np.array([0.1,-0.01,0.15,-0.05,-0.1,0.15])
    self.e_n = np.array([-0.05,0.025,-0.1,0.2,0.11,0.13])

    self.T_np1 = self.T
    self.T_n = 90.0

    self.t_np1 = 1.0
    self.t_n = 0.0

    self.u_n = 0.0
    self.p_n = 0.0
  
    # This is a rather boring baseline history state to probe, but I can't
    # think of a better way to get a "generic" history from a generic model
    self.hist_n = np.array([self.d_n] + list(self.bmodel.init_store()))
    self.x_trial = np.array([50,-25,150,-150,190,100.0] + [0.41])

    self.nsteps = 10
    self.etarget = np.array([0.1,-0.025,0.02,0.015,-0.02,-0.05])
    self.ttarget = 10.0

  def test_function(self):
    f_model = self.model.f(self.stress, self.d_np1, self.T)
    f_calcd = (self.d_np1 + self.k0) ** self.a * self.effective(self.stress) / self.W0

    self.assertTrue(np.isclose(f_model, f_calcd))


class TestCombinedDamage(unittest.TestCase, CommonScalarDamageModel,
    CommonDamagedModel):
  def setUp(self):
    self.E = 92000.0
    self.nu = 0.3

    self.s0 = 180.0
    self.Kp = 1000.0
    self.H = 1000.0

    self.elastic = elasticity.IsotropicLinearElasticModel(self.E,
        "youngs", self.nu, "poissons")

    surface = surfaces.IsoKinJ2()
    iso = hardening.LinearIsotropicHardeningRule(self.s0, self.Kp)
    kin = hardening.LinearKinematicHardeningRule(self.H)
    hrule = hardening.CombinedHardeningRule(iso, kin)

    flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

    self.bmodel = models.SmallStrainRateIndependentPlasticity(self.elastic, 
        flow)

    self.W0 = 10.0
    self.k0 = 0.0001
    self.a0 = 2.0

    self.model1 = damage.NEMLExponentialWorkDamagedModel_sd(
        self.elastic, self.W0, self.k0, 
        self.a0, self.bmodel)

    self.W02 = 10.0
    self.k02 = 0.001
    self.a02 = 1.5

    self.model2 = damage.NEMLExponentialWorkDamagedModel_sd(
        self.elastic, self.W02, self.k02, 
        self.a02, self.bmodel)

    self.model = damage.CombinedDamageModel_sd(self.elastic, 
        [self.model1, self.model2], self.bmodel)

    self.stress = np.array([100,-50.0,300.0,-99,50.0,125.0])
    self.T = 100.0
    self.d = 0.45

    self.s_np1 = self.stress
    self.s_n = np.array([-25,150,250,-25,-100,25])

    self.d_np1 = 0.5
    self.d_n = 0.4

    self.e_np1 = np.array([0.1,-0.01,0.15,-0.05,-0.1,0.15])
    self.e_n = np.array([-0.05,0.025,-0.1,0.2,0.11,0.13])

    self.T_np1 = self.T
    self.T_n = 90.0

    self.t_np1 = 1.0
    self.t_n = 0.0

    self.u_n = 0.0
    self.p_n = 0.0
  
    # This is a rather boring baseline history state to probe, but I can't
    # think of a better way to get a "generic" history from a generic model
    self.hist_n = np.array([self.d_n] + list(self.bmodel.init_store()))
    self.x_trial = np.array([50,-25,150,-150,190,100.0] + [0.41])

    self.nsteps = 10
    self.etarget = np.array([0.1,-0.025,0.02,0.015,-0.02,-0.05])
    self.ttarget = 10.0
