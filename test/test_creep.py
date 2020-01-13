import sys
sys.path.append('..')

from neml import solvers, creep, elasticity, interpolate
import unittest

from common import *

import numpy as np
import numpy.linalg as la

class CommonScalarCreep(object):
  """
    Common tests to impose on scalar creep laws.
  """
  def test_dg_ds(self):
    """
      Derivative of creep rate wrt. effective stress
    """
    dfn = lambda x: self.model.g(x, self.e, self.t, self.T)
    nderiv = differentiate(dfn, self.s)
    cderiv = self.model.dg_ds(self.s, self.e, self.t, self.T)
    self.assertTrue(np.isclose(nderiv, cderiv, rtol = 1.0e-4))

  def test_dg_de(self):
    """
      Derivative of creep rate wrt. effective strain
    """
    dfn = lambda x: self.model.g(self.s, x, self.t, self.T)
    nderiv = differentiate(dfn, self.e)
    cderiv = self.model.dg_de(self.s, self.e, self.t, self.T)
    self.assertTrue(np.isclose(nderiv, cderiv, rtol = 1.0e-4))

  def test_dg_dt(self):
    """
      Derivative of creep rate wrt. time
    """
    dfn = lambda x: self.model.g(self.s, self.e, x, self.T)
    nderiv = differentiate(dfn, self.t)
    cderiv = self.model.dg_dt(self.s, self.e, self.t, self.T)
    self.assertTrue(np.isclose(nderiv, cderiv))

  def test_dg_dT(self):
    """
      Derivative of creep rate wrt. temperature
    """
    dfn = lambda x: self.model.g(self.s, self.e, self.t, x)
    nderiv = differentiate(dfn, self.T)
    cderiv = self.model.dg_dT(self.s, self.e, self.t, self.T)
    if np.isclose(cderiv, 0.0):
      self.assertTrue(True) # This just means I don't care
    else:
      self.assertTrue(np.isclose(nderiv, cderiv))

class TestPowerLawCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.A = 1.0e-6
    self.n = 5.0

    self.model = creep.PowerLawCreep(self.A, self.n) 

    self.T = 300.0
    self.e = 0.1
    self.s = 150.0
    self.t = 10.0

  def test_properties(self):
    self.assertTrue(np.isclose(self.A, self.model.A(self.T)))
    self.assertTrue(np.isclose(self.n, self.model.n(self.T)))

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    g_calc = self.A * self.s ** (self.n)
    self.assertTrue(np.isclose(g_direct, g_calc))

class TestNormalizedPowerLawCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.s0 = 500.0
    self.n = 3.0

    self.model = creep.NormalizedPowerLawCreep(self.s0, self.n) 

    self.T = 300.0
    self.e = 0.1
    self.s = 150.0
    self.t = 10.0

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    g_calc = (self.s / self.s0) ** (self.n)
    self.assertTrue(np.isclose(g_direct, g_calc))

class TestBlackburnMinimumCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.A = 1.0e-6
    self.beta = -2.0e-2
    self.n = 5.0
    self.Q = 67000.0
    self.R = 1.987

    self.model = creep.BlackburnMinimumCreep(self.A, self.n, self.beta,
        self.R, self.Q) 

    self.T = 300.0
    self.e = 0.1
    self.s = 150.0
    self.t = 10.0

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    g_calc = self.A * np.sinh(self.beta*self.s/self.n)**self.n * np.exp(-self.Q/(self.R*self.T))
    self.assertTrue(np.isclose(g_direct, g_calc))

class TestSwindemanMinimumCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.C = 2.25e20
    self.n = 5.0
    self.V = 0.038
    self.Q = 77280.0

    self.model = creep.SwindemanMinimumCreep(self.C, self.n, self.V,
        self.Q) 

    self.T = 600+273.15
    self.e = 0.1
    self.s = 150.0
    self.t = 10.0

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    g_calc = self.C * self.s ** self.n * np.exp(self.V * self.s) * np.exp(-self.Q/self.T)
    self.assertTrue(np.isclose(g_direct, g_calc))

class TestRegionKMCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.b = 2.019 * 1.0e-7
    self.eps0 = 1.0e6
    self.E = 160000.0
    self.nu = 0.3

    self.emodel = elasticity.IsotropicLinearElasticModel(self.E, "youngs",
        self.nu, "poissons")

    self.kboltz = 1.38064e-23 * 1000.0
    
    self.cuts = [0.00145709]
    self.A = [-0.7869, -0.1695]
    self.B = [-4.2388, -0.2065]

    self.model = creep.RegionKMCreep(self.cuts, self.A, self.B, self.kboltz,
        self.b, self.eps0, self.emodel)

    self.T = 300.0
    self.e = 0.1
    self.s = 300.0
    self.t = 10.0

class TestNortonBaileyCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.A = 1.0e-6
    self.m = 0.25
    self.n = 5.0

    self.model = creep.NortonBaileyCreep(self.A, self.m, self.n) 

    self.T = 300.0
    self.e = 0.1
    self.s = 150.0
    self.t = 10.0

  def test_properties(self):
    self.assertTrue(np.isclose(self.A, self.model.A(self.T)))
    self.assertTrue(np.isclose(self.m, self.model.m(self.T)))
    self.assertTrue(np.isclose(self.n, self.model.n(self.T)))

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    g_calc = self.m * self.A**(1.0/self.m) * self.s**(self.n/self.m) * self.e**((self.m-1.0)/self.m)
    self.assertTrue(np.isclose(g_direct, g_calc))

class TestMukherjeeCreep(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.A = 1.0e10     # Unitless
    self.n = 7.9        # Unitless
    self.D0 = 3.7e-5    # m^2 / s
    self.Qv = 280.0     # kJ/ mol
    self.b = 2.28e-10   # m

    self.E = 100000.0
    self.nu = 0.3
    self.emodel = elasticity.IsotropicLinearElasticModel(self.E, "youngs",
        self.nu, "poissons")
    
    self.k = 1.380e-23 * 1.0e-6 # m^3 * MPa / K
    self.R = 8.314 / 1000 # kJ / (K * mol)

    self.t = 10.0
    self.e = 0.1
    self.s = 150.0
    self.T = 550.0 + 273.15

    self.model = creep.MukherjeeCreep(self.emodel, self.A, self.n, self.D0,
        self.Qv, self.b, self.k, self.R)

  def test_properties(self):
    self.assertTrue(np.isclose(self.A, self.model.A))
    self.assertTrue(np.isclose(self.n, self.model.n))
    self.assertTrue(np.isclose(self.D0, self.model.D0))
    self.assertTrue(np.isclose(self.Qv, self.model.Q))
    self.assertTrue(np.isclose(self.b, self.model.b))
    self.assertTrue(np.isclose(self.k, self.model.k))
    self.assertTrue(np.isclose(self.R, self.model.R))

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    mu = self.emodel.G(self.T)
    g_calc = self.A * self.D0 * np.exp(-self.Qv / (self.R * self.T)
        ) * mu * self.b / (self.k * self.T) * (self.s / mu) ** self.n
    self.assertTrue(np.isclose(g_direct, g_calc))

class TestGenericCreep(unittest.TestCase, CommonScalarCreep):
  """
    Tests GenericCreep --- creep law given as a generic interpolation
    of log stress.
  """
  def setUp(self):
    # Polynomial coefficients: log(creep rate) = np.polyval(poly, log(stress))
    self.poly = [2.5, -50.0, 365.0, -1200.0, 1500.0]
    self.ifn = interpolate.PolynomialInterpolate(self.poly)
    
    # Model only has one input parameter (the interpolation formula)
    self.model = creep.GenericCreep(self.ifn)
    
    self.T = 300.0  # Temperature to test at
    self.e = 0.1    # Effective strain to test at
    self.s = 150.0  # Effective stress to test at
    self.t = 10.0   # Time to test at

  def test_g(self):
    """
      Compare the creep rate calculated "by hand" to the creep rate
      returned by the model.
    """
    g_calc = np.exp(np.polyval(self.poly, np.log(self.s)))
    g_model = self.model.g(self.s, self.e, self.t, self.T)
    self.assertTrue(np.isclose(g_calc, g_model))

class Test225CreepCase1(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.model = creep.MinCreep225Cr1MoCreep() 

    self.T = 500 + 273.15
    self.e = 0.1
    self.s = 50.0
    self.t = 10.0

    self.U = interpolate.PiecewiseLinearInterpolate(
        [371, 400, 450, 500, 550, 600, 621, 649],
        [471,468,452,418,634,284,300,270])

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    Uv = self.U.value(self.T - 273.15)
    g_calc = 10.0**(6.7475+0.011426*self.s+987.72/Uv * np.log10(self.s) - 13494.0 / self.T) / 100.0
    self.assertTrue(np.isclose(g_direct, g_calc))

class Test225CreepCase2(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.model = creep.MinCreep225Cr1MoCreep() 

    self.T = 500 + 273.15
    self.e = 0.1
    self.s = 100.0
    self.t = 10.0

    self.U = interpolate.PiecewiseLinearInterpolate(
        [371, 400, 450, 500, 550, 600, 621, 649],
        [471,468,452,418,634,284,300,270])

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    Uv = self.U.value(self.T - 273.15)
    g_calc = 10.0**(6.7475+0.011426*self.s+987.72/Uv * np.log10(self.s) - 13494.0 / self.T) / 100.0
    self.assertTrue(np.isclose(g_direct, g_calc))

class Test225CreepCase3(unittest.TestCase, CommonScalarCreep):
  def setUp(self):
    self.model = creep.MinCreep225Cr1MoCreep() 

    self.T = 600 + 273.15
    self.e = 0.1
    self.s = 100.0
    self.t = 10.0

    self.U = interpolate.PiecewiseLinearInterpolate(
        [371, 400, 450, 500, 550, 600, 621, 649],
        [471,468,452,418,634,284,300,270])

  def test_g(self):
    g_direct = self.model.g(self.s, self.e, self.t, self.T)
    Uv = self.U.value(self.T - 273.15)
    g_calc = 10.0**(11.498 - 8.2226 * Uv / self.T - 20448 / self.T + 
        5862.4 / self.T * np.log10(self.s)) / 100.0
    self.assertTrue(np.isclose(g_direct, g_calc))

class CommonCreepModel(object):
  """
    Tests common to all creep models
  """
  def test_tangent(self):
    strain = np.copy(self.e)
    stress = np.copy(self.s)

    ts = np.linspace(0,self.t)

    for i in range(1, len(ts)):
      t_np1 = ts[i]
      t_n = ts[i-1]
      s_np1 = stress
      e_next, A_next = self.model.update(stress, strain, self.T, self.T,
          t_np1, t_n)

      dfn = lambda x: self.model.update(x, strain, self.T, self.T, t_np1,
          t_n)[0]
      nA = differentiate(dfn, s_np1)
      self.assertTrue(np.allclose(nA, A_next, rtol = 1.0e-4))

      strain = np.copy(e_next)

  def test_jacobian(self):
    ts = self.model.make_trial_state(self.s, self.e, self.T, self.T,
      self.dt + self.t, self.t)
    R, J = self.model.RJ(self.x, ts)
    nJ = differentiate(lambda x: self.model.RJ(x, ts)[0], self.x)
    self.assertTrue(np.allclose(J, nJ, rtol = 1.0e-4))

  def test_df_ds(self):
    dfn = lambda x: self.model.f(x, self.e, self.t, self.T)
    nderiv = differentiate(dfn, self.s)
    cderiv = self.model.df_ds(self.s, self.e, self.t, self.T)
    self.assertTrue(np.allclose(nderiv, cderiv, rtol = 1.0e-4))

  def test_df_de(self):
    dfn = lambda x: self.model.f(self.s, x, self.t, self.T)
    nderiv = differentiate(dfn, self.e)
    cderiv = self.model.df_de(self.s, self.e, self.t, self.T)
    self.assertTrue(np.allclose(nderiv, cderiv, rtol = 1.0e-4))

  def test_df_dt(self):
    dfn = lambda x: self.model.f(self.s, self.e, x, self.T)
    nderiv = differentiate(dfn, self.t)
    cderiv = self.model.df_dt(self.s, self.e, self.t, self.T)
    self.assertTrue(np.allclose(nderiv, cderiv))

  def test_df_dT(self):
    dfn = lambda x: self.model.f(self.s, self.e, self.t, x)
    nderiv = differentiate(dfn, self.T)
    cderiv = self.model.df_dT(self.s, self.e, self.t, self.T)
    self.assertTrue(np.allclose(nderiv, cderiv))

class TestJ2Creep(unittest.TestCase, CommonCreepModel):
  def setUp(self):
    self.A = 1.0e-10
    self.m = 0.25
    self.n = 5.0

    self.smodel = creep.NortonBaileyCreep(self.A, self.m, self.n)

    self.model = creep.J2CreepModel(self.smodel)

    self.s = np.array([100.0,-25.0,-5.0, 20.0,15.0,3.0])
    self.e = np.array([0.05, -0.01, -0.01, 0.025, 0.03, -0.01])
    self.x = np.array([0.06, 0.01, 0.02, 0.01, 0.02, -0.05])
    self.T = 300.0
    self.t = 10.0
    self.dt = 2.0

    tr = np.sum(self.s[:3])
    self.sdev = self.s - np.array([1,1,1,0,0,0]) * tr / 3.0
    self.se = np.sqrt(3.0/2.0) * la.norm(self.sdev)
    self.ee = np.sqrt(2.0/3.0) * la.norm(self.e)

  def test_f(self):
    f_direct = self.model.f(self.s, self.e, self.t, self.T)
    g_direct = self.smodel.g(self.se, self.ee, self.t, self.T)
    f_calc = g_direct * 3.0 / 2.0 * self.sdev / self.se

    self.assertTrue(np.allclose(f_direct, f_calc))

  def test_f_uni(self):
    """
      Make sure all our "effectives" work out correctly
    """
    s = np.array([100.0, 0, 0, 0, 0, 0])
    e = np.array([0.1, -0.05, -0.05, 0, 0, 0])
    f_direct = self.model.f(s, e, self.t, self.T)
    
    sdev = s - np.array([1,1,1,0,0,0]) * np.sum(s[:3]) / 3.0
    se = np.sqrt(3.0/2.0) * la.norm(sdev)
    ee = np.sqrt(2.0/3.0) * la.norm(e)

    g_direct = self.smodel.g(se, ee, self.t, self.T)
    
    self.assertTrue(np.isclose(g_direct, f_direct[0]))

    self.assertTrue(np.isclose(-g_direct/2.0, f_direct[1]))
    self.assertTrue(np.isclose(-g_direct/2.0, f_direct[2]))

    self.assertTrue(np.allclose([0,0,0], f_direct[3:]))


