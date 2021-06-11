from neml import models, interpolate, elasticity, history
from neml.cp import hucocks
from neml.math import rotations, tensors

import unittest
import numpy as np
import scipy.interpolate as inter

from common import differentiate_new

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
    print(c)
    print(self.kboltz, T, self.vm, np.prod(c), np.prod(self.ceq(T)))
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
