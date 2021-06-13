from neml import models, interpolate, elasticity, history
from neml.cp import hucocks, crystallography, sliprules
from neml.math import rotations, tensors

import unittest
import numpy as np
import scipy.interpolate as inter

from common import differentiate_new
from test_slipharden import CommonSlipHardening 

class PrecipitationModel:
  """
    Implementation that can/will mirror a NEML implementation
  """
  def __init__(self, c0, cp, ceq, am, N0, Vm, chi, D0, Q0, 
      Cf, kboltz = 1.3806485e-23, R = 8.31462, Na = 6.02e23,
      use = 0):
    """
      Parameters:
        c0:     initial matrix concentration, as a function of temperature
        cp:     precipitate concentration, as a function of temperature
        ceq:    equilibrium matrix concentration, as a function of temperature
        am:     lattice parameter
        No:     potential nucleation sites per volume
        Vm:     species molar volume
        chi:    surface energy
        D0:     diffusion prefactor
        Q0:     diffusion activation energy
        Cf:     coarsening reduction factor

      Additional Parameters:
        kboltz:     boltzman constant in the right units
        R:          gas constant in the right units
        Na:         Avagadro's number
        use:        which component of the concentration to use to drive kinetics
    """
    self.c0 = c0
    self.cp = cp
    self.ceq = ceq
    self.am = am
    self.N0 = N0
    self.Vm = Vm
    # Molecular volume
    self.vm = self.Vm / Na
    self.chi = chi
    self.D0 = D0
    self.Q0 = Q0
    
    self.Cf = Cf

    self.kboltz = kboltz
    self. R = R

    self.use = use

  def __call__(self, t, y, T):
    """
      The actual ODE (the h_dot in NEML)

      We can get away with 2 variables, which is nice

      Stack the vector as [r, N]

      Parameters:
        t:      time...
        y:      state...
        T:      temperature...
    """
    ydot = np.zeros((3,))
    ydot[0] = self.f_dot(y[0], y[1], y[2], T)
    ydot[1] = self.r_dot(y[0], y[1], y[2], T)
    ydot[2] = self.N_dot(y[0], y[1], y[2], T)

    return ydot

  def D(self, T):
    """
      Full temperature-dependent diffusivity
    """
    return self.D0 * np.exp(-self.Q0 / (self.R * T))

  def Gv(self, c, T):
    """
      Gibbs free energy driving precipitation

      Parameters:
        c:      vector of current concentrations
        T:      temperature
    """
    return -self.kboltz * T / self.vm * np.log(np.prod(c) / np.prod(self.ceq(T)))

  def c(self, f, T):
    """
      Calculate the current concentration vector given the current
      volume fraction

      Parameters:
        f:      volume fraction
        T:      temperature
    """
    return (self.c0(T) - f * self.cp(T)) / (1.0 - f)

  def f_dot(self, f, r, N, T):
    return 4.0*np.pi/3 * (self.N_dot(f, r, N, T) * r**3.0 + N * 3.0 * r**2.0
        * self.r_dot(f, r, N, T))

  def r_dot(self, f, r, N, T):
    """
      Radius growth rate

      Parameters:
        r:      current radius
        N:      current number density
        T:      current temperature
    """
    c = self.c(f, T)
    D = self.D(T)

    # Branch 1: nucleation
    if np.all(c >= self.ceq(T)):
      Gv = self.Gv(c, T)
      rc = -2.0*self.chi / Gv

      rt = D/r * (c[self.use] - self.ceq(T)[self.use]) / (self.cp(T)[self.use] - self.ceq(T)[self.use]
          ) + self.N_dot(f, r, N, T) / N * (rc - r)
      return rt
    # Branch 2: ripening
    else:
      K = self.Cf(T) * 8.0 * self.chi * self.Vm * D * c[self.use] / (
          9.0 * self.R * T)
      rt = K / (3.0*r**2.0)
      return rt

  def N_dot(self, f, r, N, T):
    """
      Number density growth rate

      Parameters:
        r:      current radius
        N:      current number density
        T:      current temperature
    """
    c = self.c(f, T)

    # Branch 1: nucleation
    if np.all(c >= self.ceq(T)):
      Gv = self.Gv(c, T)
      Gstar = 16*np.pi*self.chi**3.0/(3.0*Gv**2.0)
      D = self.D(T)
      ZB = 2*self.vm * D  * c[self.use]/ self.am**4.0 * np.sqrt(
          self.chi/(self.kboltz * T))
      Nt = self.N0 * ZB * np.exp(-Gstar / (self.kboltz * T))
      return Nt
    # Branch 2: ripening
    else:
      return -3.0 * N / r * self.r_dot(f, r, N, T)

  def ic(self, T):
    """
      Initial condition for both variables
    """
    return np.array([0,1.0e-20,1.0e-6])


class TestPercipitationModel(unittest.TestCase):
  def setUp(self):
    # Data from paper
    Ts = np.array([500.0,550.0,600.0,650.0]) + 273.15

    # Setup for [Cr,C] <-> M23C6
    self.am = 3.6e-10
    self.N0 = 1.0e13
    self.Vm = 6e-6
    self.chi = 0.3
    self.D0 = 1.5e-4
    self.Q0 = 240e3
    self.c0_py = lambda T: np.array([16.25,0.0375])/100.0
    self.cp_py = lambda T: inter.interp1d(
        Ts,
        np.array([[69.85,5.13],[69.05,5.13],[68.32,5.13],[67.52,5.13]])/100.0, axis = 0,
        fill_value = "extrapolate")(T)
    self.ceq_py = lambda T: inter.interp1d(
        Ts,
        np.array([[15.64,7.25e-6],[15.69,2.92e-5],[15.75,9.48e-5],[15.83,2.97e-4]])/100.0, axis = 0,
        fill_value = "extrapolate")(T)
    self.Cf_py = lambda T: inter.interp1d(Ts, [1.0,1.0,0.3,0.03],
        fill_value = "extrapolate")(T)


    self.model_py = PrecipitationModel(self.c0_py, self.cp_py, self.ceq_py,
        self.am, self.N0, self.Vm, self.chi, self.D0, self.Q0, self.Cf_py)

    self.c0_neml = [interpolate.ConstantInterpolate(16.25/100.0),
        interpolate.ConstantInterpolate(0.0375/100.0)]
    self.cp_neml = [interpolate.PiecewiseLinearInterpolate(list(Ts),
      [69.85/100, 69.05/100, 68.32/100, 67.52/100]),
      interpolate.PiecewiseLinearInterpolate(list(Ts),
        [5.13/100, 5.13/100, 5.13/100, 5.13/100])]
    self.ceq_neml = [interpolate.PiecewiseLinearInterpolate(list(Ts),
      [15.64/100,15.69/100,15.75,100,15.83/100]),
      interpolate.PiecewiseLinearInterpolate(list(Ts),
        [7.25e-6/100, 2.92e-5/100, 9.48e-5/100, 2.97e-4/100])]
    self.Cf_neml = interpolate.PiecewiseLinearInterpolate(list(Ts), 
        [1.0, 1.0, 0.3, 0.03])

    self.model_neml = hucocks.HuCocksPrecipitationModel(self.c0_neml, self.cp_neml, self.ceq_neml, 
        self.am, self.N0, self.Vm, self.chi, self.D0, self.Q0, self.Cf_neml)
   
  def make_state(self, f, r, N):
    hist = history.History()

    hist.add_scalar("f")
    hist.set_scalar("f", f)
    hist.add_scalar("r")
    hist.set_scalar("r", r)
    hist.add_scalar("N")
    hist.set_scalar("N", N)

    return hist

  def test_fdot_nucleation(self):
    """
      Test f_dot in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    py_v = self.model_py.f_dot(f, r, N, T)
    ne_v = self.model_neml.f_rate(self.make_state(f, r, N), T)
    self.assertAlmostEqual(py_v, ne_v)

  def test_df_df_nucleation(self):
    """
      Test df_df in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.f_rate(self.make_state(x, r, N), T), f)
    exact = self.model_neml.df_df(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_df_dr_nucleation(self):
    """
      Test df_dr in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.f_rate(self.make_state(f, x, N), T), r)
    exact = self.model_neml.df_dr(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_df_dN_nucleation(self):
    """
      Test df_dN in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.f_rate(self.make_state(f, r, x), T), N)
    exact = self.model_neml.df_dN(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_fdot_ostwald(self):
    """
      Test f_dot in the ripening regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    py_v = self.model_py.f_dot(f, r, N, T)
    ne_v = self.model_neml.f_rate(self.make_state(f, r, N), T)
    self.assertAlmostEqual(py_v, ne_v)

  def test_df_df_ostwald(self):
    """
      Test df_df in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.f_rate(self.make_state(x, r, N), T), f)
    exact = self.model_neml.df_df(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_df_dr_ostwald(self):
    """
      Test df_dr in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.f_rate(self.make_state(f, x, N), T), r)
    exact = self.model_neml.df_dr(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_df_dN_ostwald(self):
    """
      Test df_dN in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.f_rate(self.make_state(f, r, x), T), N)
    exact = self.model_neml.df_dN(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_rdot_nucleation(self):
    """
      Test r_dot in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    py_v = self.model_py.r_dot(f, r, N, T)
    ne_v = self.model_neml.r_rate(self.make_state(f, r, N), T)
    self.assertAlmostEqual(py_v, ne_v)

  def test_dr_df_nucleation(self):
    """
      Test dr_df in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.r_rate(self.make_state(x, r, N), T), f)
    exact = self.model_neml.dr_df(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_dr_dr_nucleation(self):
    """
      Test dr_dr in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.r_rate(self.make_state(f, x, N), T), r)
    exact = self.model_neml.dr_dr(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_dr_dN_nucleation(self):
    """
      Test dr_dN in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.r_rate(self.make_state(f, r, x), T), N)
    exact = self.model_neml.dr_dN(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_rdot_ostwald(self):
    """
      Test r_dot in the ripening regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    py_v = self.model_py.r_dot(f, r, N, T)
    ne_v = self.model_neml.r_rate(self.make_state(f, r, N), T)
    self.assertAlmostEqual(py_v, ne_v)

  def test_dr_df_ostwald(self):
    """
      Test dr_df in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.r_rate(self.make_state(x, r, N), T), f)
    exact = self.model_neml.dr_df(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_dr_dr_ostwald(self):
    """
      Test dr_dr in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.r_rate(self.make_state(f, x, N), T), r)
    exact = self.model_neml.dr_dr(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_dr_dN_ostwald(self):
    """
      Test dr_dN in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.r_rate(self.make_state(f, r, x), T), N)
    exact = self.model_neml.dr_dN(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_Ndot_nucleation(self):
    """
      Test N_dot in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    py_v = self.model_py.N_dot(f, r, N, T)
    ne_v = self.model_neml.N_rate(self.make_state(f, r, N), T)
    self.assertAlmostEqual(py_v, ne_v)

  def test_dN_df_nucleation(self):
    """
      Test dN_df in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.N_rate(self.make_state(x, r, N), T), f)
    exact = self.model_neml.dN_df(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_dN_dr_nucleation(self):
    """
      Test dN_dr in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.N_rate(self.make_state(f, x, N), T), r)
    exact = self.model_neml.dN_dr(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_dN_dN_nucleation(self):
    """
      Test dN_dN in the nucleation regime
    """
    T = 550.0 + 273.15
    f = 0.0035
    r = 34e-9
    N = f / (4.0/3.0*np.pi*r**3.0)

    num = differentiate_new(lambda x: self.model_neml.N_rate(self.make_state(f, r, x), T), N)
    exact = self.model_neml.dN_dN(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_Ndot_ostwald(self):
    """
      Test N_dot in the ripening regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    py_v = self.model_py.N_dot(f, r, N, T)
    ne_v = self.model_neml.N_rate(self.make_state(f, r, N), T)
    self.assertAlmostEqual(py_v, ne_v)

  def test_dN_df_ostwald(self):
    """
      Test dN_df in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.N_rate(self.make_state(x, r, N), T), f)
    exact = self.model_neml.dN_df(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_dN_dr_ostwald(self):
    """
      Test dN_dr in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.N_rate(self.make_state(f, x, N), T), r)
    exact = self.model_neml.dN_dr(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-4))

  def test_dN_dN_ostwald(self):
    """
      Test dN_dN in the ostwald regime
    """
    T = 550.0 + 273.15
    f = 0.073
    r = 95e-9
    N = f / (4.0/3.0*np.pi*r**3.0) 

    num = differentiate_new(lambda x: self.model_neml.N_rate(self.make_state(f, r, x), T), N)
    exact = self.model_neml.dN_dN(self.make_state(f, r, N), T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

  def test_c(self):
    """
      Test the conversion from f to c
    """
    T = 550.0 + 273.15
    f = 0.0035

    py_c = self.model_py.c(f, T)
    ne_c = self.model_neml.c(f, T)

    self.assertTrue(np.allclose(py_c, ne_c))

  def test_dc_df(self):
    """
      Test the derivative of the conversion from f to c
    """
    T = 550.0 + 273.15
    f = 0.0035

    num = differentiate_new(lambda x: np.array(self.model_neml.c(x, T)), np.array([f])).flatten()
    exact = self.model_neml.dc_df(f, T)

    self.assertTrue(np.allclose(num, exact))

  def test_Gv(self):
    """
      Test the calculation of the Gibbs free energy driving precipitation
    """
    T = 550.0 + 273.15
    f = 0.0035

    py_c = self.model_py.Gv(self.model_py.c(f, T), T)
    ne_c = self.model_neml.Gv(f, T)

    self.assertAlmostEqual(py_c, ne_c)

  def test_dG_df(self):
    """
      Test the derivative of the Gibbs free energy driving precipitation
    """
    T = 550.0 + 273.15
    f = 0.0035

    num = differentiate_new(lambda x: self.model_neml.Gv(x, T), f)
    exact = self.model_neml.dG_df(f, T)

    print(num, exact)
    self.assertTrue(np.isclose(num, exact, rtol=1e-5))

class TestDislocationSpacingHardening(unittest.TestCase, CommonSlipHardening):
  def setUp(self):
    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))
    
    self.nslip = self.L.ntotal

    self.current = 1.15e-7

    self.H = history.History()
    for i in range(self.nslip):
      self.H.add_scalar("spacing_"+str(i))
      self.H.set_scalar("spacing_"+str(i), self.current)

    self.T = 600 + 273.15
    
    self.J1 = 2e14
    self.J2 = 3.3e14
    self.K = 2.56e-30
    self.L0 = 3.1623e-7
    self.b = 2.5e-10
    self.a = 0.35
    self.G = 57633.58
    
    self.model = hucocks.DislocationSpacingHardening(self.J1, self.J2, self.K, self.L0, self.a, 
        self.b, self.G, self.L)

    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)

    self.fixed = history.History()

  def test_hist_to_tau(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        model = self.model.hist_to_tau(g, i, self.H, self.L, self.T,
            self.fixed)
        should = self.a * self.G * self.b / self.current
        self.assertAlmostEqual(model, should)


  def test_definition(self):
    hrate = self.model.hist(self.S, self.Q, self.H, self.L, self.T, self.sliprule,
        self.fixed)
    should = np.zeros((self.nslip,))

    for i in range(self.nslip):
      for g in range(self.L.ngroup):
        for k in range(self.L.nslip(g)):
          j = self.L.flat(g,k)
          if i == j:
            c = self.J1
          else:
            c = self.J2
          should[i] -= self.current**3.0 * c * np.abs(
              self.sliprule.slip(g, k, self.S, self.Q, self.H, self.L, self.T, self.fixed))
      should[i] += self.K / self.current**3.0
    
    self.assertTrue(np.allclose(hrate, should))

class TestHuCocksHardening(unittest.TestCase, CommonSlipHardening):
  def setUp(self):
    Ts = np.array([500.0,550.0,600.0,650.0]) + 273.15

    self.L = crystallography.CubicLattice(1.0)
    self.L.add_slip_system([1,1,0],[1,1,1])
    
    self.Q = rotations.Orientation(35.0,17.0,14.0, angle_type = "degrees")
    self.S = tensors.Symmetric(np.array([
      [100.0,-25.0,10.0],
      [-25.0,-17.0,15.0],
      [10.0,  15.0,35.0]]))

    self.T = 600 + 273.15
    
    self.J1 = 2e14
    self.J2 = 3.3e14
    self.K = 2.56e-30
    self.L0 = 3.1623e-7
    self.b = 2.5e-10
    self.ad = 0.35
    self.G = interpolate.PiecewiseLinearInterpolate(
        list(Ts),
        [61068, 59541.0, 57633.6, 55725.2])
    
    self.dmodel = hucocks.DislocationSpacingHardening(self.J1, self.J2, self.K, self.L0, self.ad, 
        self.b, self.G, self.L)

    # Setup for [Cr,C] <-> M23C6
    self.am_car = 3.6e-10
    self.N0_car = 1.0e13
    self.Vm_car = 6e-6
    self.chi_car = 0.3
    self.D0_car = 1.5e-4
    self.Q0_car = 240e3
    self.c0_car = [interpolate.ConstantInterpolate(16.25/100.0),
        interpolate.ConstantInterpolate(0.0375/100.0)]
    self.cp_car = [interpolate.PiecewiseLinearInterpolate(list(Ts),
      [69.85/100, 69.05/100, 68.32/100, 67.52/100]),
      interpolate.PiecewiseLinearInterpolate(list(Ts),
        [5.13/100, 5.13/100, 5.13/100, 5.13/100])]
    self.ceq_car = [interpolate.PiecewiseLinearInterpolate(list(Ts),
      [15.64/100,15.69/100,15.75,100,15.83/100]),
      interpolate.PiecewiseLinearInterpolate(list(Ts),
        [7.25e-6/100, 2.92e-5/100, 9.48e-5/100, 2.97e-4/100])]
    self.Cf_car = interpolate.PiecewiseLinearInterpolate(list(Ts), 
        [1.0, 1.0, 0.3, 0.03])

    self.carbide = hucocks.HuCocksPrecipitationModel(self.c0_car, self.cp_car, self.ceq_car, 
        self.am_car, self.N0_car, self.Vm_car, self.chi_car, self.D0_car,
        self.Q0_car, self.Cf_car) 

    self.am_laves = 3.6e-10
    self.N0_laves = 5e14
    self.Vm_laves = 2e-6
    self.chi_laves = 0.25
    self.D0_laves = 7.4e-4
    self.Q0_laves = 283e3
    self.c0_laves = [2.33/100.0]
    self.cp_laves = [50.0/100.0]
    self.ceq_laves = [interpolate.PiecewiseLinearInterpolate(list(Ts),
      [0.25/100,0.46/100.0,0.76/100.0,1.16/100.0])]
    self.Cf_laves = 1.0

    self.laves = hucocks.HuCocksPrecipitationModel(self.c0_laves, self.cp_laves, self.ceq_laves, 
        self.am_laves, self.N0_laves, self.Vm_laves, self.chi_laves, self.D0_laves,
        self.Q0_laves, self.Cf_laves) 

    self.ap = 0.84
    self.ac = 0.000457

    self.model = hucocks.HuCocksHardening(self.dmodel, [self.carbide, self.laves], self.ap, self.ac, self.b,
        self.G)

    self.g0 = 1.0
    self.n = 3.0
    self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)
    
    self.nslip = self.L.ntotal

    self.current = 1.15e-7
    self.f_carbide = 0.073
    self.r_carbide = 95e-9
    self.N_carbide = self.f_carbide / (4.0/3.0*np.pi*self.r_carbide**3.0)
    self.f_laves = 0.00962
    self.r_laves = 41.2e-9
    self.N_laves = self.f_laves / (4.0/3.0*np.pi*self.r_laves**3.0) 

    self.H = history.History()
    for i in range(self.nslip):
      self.H.add_scalar("spacing_"+str(i))
      self.H.set_scalar("spacing_"+str(i), self.current)
    
    self.H.add_scalar("f_0")
    self.H.set_scalar("f_0", self.f_carbide)
    self.H.add_scalar("r_0")
    self.H.set_scalar("r_0", self.r_carbide)
    self.H.add_scalar("N_0")
    self.H.set_scalar("N_0", self.N_carbide)

    self.H.add_scalar("f_1")
    self.H.set_scalar("f_1", self.f_laves)
    self.H.add_scalar("r_1")
    self.H.set_scalar("r_1", self.r_laves)
    self.H.add_scalar("N_1")
    self.H.set_scalar("N_1", self.N_laves)

    self.fixed = history.History()

  def test_history(self):
    """
      Paranoid check on the history state, given this is a kinda complicated model
    """
    H1 = history.History()
    self.model.populate_history(H1)
    self.model.init_history(H1)

    self.assertEqual(len(np.array(H1)), 12+3+3)
    
    should = np.array([self.L0]*12+[0, 1.0e-20, 1.0e-6]*2) 
    
    print(np.array(H1))
    print(should)

    self.assertTrue(np.allclose(np.array(H1), should))

  def test_hist_to_tau(self):
    for g in range(self.L.ngroup):
      for i in range(self.L.nslip(g)):
        model = self.model.hist_to_tau(g, i, self.H, self.L, self.T,
            self.fixed)

        tau_d = self.dmodel.hist_to_tau(g, i, self.H, self.L, self.T,
            self.fixed)
        Na = 2.0*self.H.get_scalar("r_0") * self.H.get_scalar("N_0") + 2.0*self.H.get_scalar("r_1") * self.H.get_scalar("N_1") 
        c = np.sum(self.carbide.c(self.H.get_scalar("f_0"), self.T))/self.carbide.vm + np.sum(
            self.laves.c(self.H.get_scalar("f_1"), self.T)) / self.laves.vm
        tau_p = self.ap * self.G.value(self.T) * self.b * np.sqrt(Na)
        tau_c = self.ac * self.G.value(self.T) * self.b * np.sqrt(c*self.b)
        should = np.sqrt(tau_d**2.0 + tau_p**2.0) + tau_c
        self.assertAlmostEqual(model, should)

  def test_definition(self):
    """
      All concatenated together
    """
    assem = np.array(list(np.array(self.dmodel.hist(self.S, self.Q, self.H, 
      self.L, self.T, self.sliprule, self.fixed))) + [
          self.carbide.f_rate(self.H, self.T),
          self.carbide.r_rate(self.H, self.T),
          self.carbide.N_rate(self.H, self.T),
          self.laves.f_rate(self.H, self.T),
          self.laves.r_rate(self.H, self.T),
          self.laves.N_rate(self.H, self.T)])
    fmodel = self.model.hist(self.S, self.Q, self.H, self.L, self.T,
        self.sliprule, self.fixed)

    self.assertTrue(np.allclose(assem, fmodel))
