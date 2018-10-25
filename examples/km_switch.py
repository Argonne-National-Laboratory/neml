#!/usr/bin/env python

import sys
sys.path.append('..')
from neml import solvers, interpolate, models, hardening, elasticity, visco_flow, general_flow, surfaces, drivers, calibrate, ri_flow

import matplotlib.pyplot as plt
import numpy as np

def rate_jump_test(model, erates, strains, T, nsteps = 50,
    sdir = np.array([1,0,0,0,0,0])):
  """
    Run a strain rate jump test.
  """
  driver = drivers.Driver_sd(model)
  strain = []
  stress = []

  for erate, de in zip(erates, strains):
    e_inc = de / nsteps
    for i in range(nsteps):
      einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T)
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

  return np.array(strain), np.array(stress)

def temp_jump_test(model, erate, strains, Ts, nsteps = 50, 
    sdir = np.array([1,0,0,0,0,0])):
  """
    Run a temperature jump test.
  """
  driver = drivers.Driver_sd(model)
  strain = []
  stress = []

  for T, de in zip(Ts, strains):
    e_inc = de / nsteps
    for i in range(nsteps):
      einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, T)
      strain.append(np.dot(driver.strain_int[-1], sdir))
      stress.append(np.dot(driver.stress_int[-1], sdir))

  return np.array(strain), np.array(stress)

if __name__ == "__main__":
  # Fully-defined perfectly plastic model
  Epoly = [-78.2759, 236951.0]
  nu = 0.3
  A = -9.6187
  B = -1.4819
  C = -5.0486
  g0 = 0.3708

  b = 0.248 * 1.0e-6
  kboltz = 1.38064e-23 * 1000.0
  eps0 = 1.0e10

  # Temperature range over which to consider (K)
  Tmin = 550.0
  Tmax = 950.0
  Trange = np.linspace(Tmin, Tmax)

  # Elastic
  E_m = interpolate.PolynomialInterpolate(Epoly)
  nu_m = interpolate.ConstantInterpolate(nu)
  elastic_m = elasticity.IsotropicLinearElasticModel(E_m, "youngs", nu_m,
      "poissons")

  # Rate sensitivity interpolates values
  mu_values = np.array([elastic_m.G(T) for T in Trange])
  n_values = -mu_values*b**3.0 / (kboltz * Trange * A)
  eta_values = np.exp(B) * eps0 ** (kboltz * Trange * A / (mu_values * b**3.0)) * mu_values

  # Rate independent interpolate values
  flow_stress = mu_values * np.exp(C)
  
  # Common objects
  surface = surfaces.IsoKinJ2()
  hmodulus = interpolate.PolynomialInterpolate([-10.0, 12000.0])
  #hmodulus = interpolate.ConstantInterpolate(1000.0)

  # Setup visco model
  n_interp = interpolate.PiecewiseLinearInterpolate(list(Trange), list(n_values))
  eta_interp = interpolate.PiecewiseLinearInterpolate(list(Trange), list(eta_values))
  eta_m = visco_flow.ConstantFluidity(eta_interp)

  iso_rd = hardening.LinearIsotropicHardeningRule(
      interpolate.ConstantInterpolate(0.0),
      hmodulus)
  hard_rd = hardening.Chaboche(iso_rd,
      [interpolate.ConstantInterpolate(0.0)], 
      [hardening.ConstantGamma(interpolate.ConstantInterpolate(0.0))],
      [interpolate.ConstantInterpolate(0.0)],
      [interpolate.ConstantInterpolate(1.0)])

  visco_rd = visco_flow.ChabocheFlowRule(surface, hard_rd, eta_m, n_interp) 
  general_rd = general_flow.TVPFlowRule(elastic_m, visco_rd)

  rate_dependent = models.GeneralIntegrator(elastic_m, general_rd)

  # Setup rate independent
  sy_interp = interpolate.PiecewiseLinearInterpolate(list(Trange), list(flow_stress))
  iso_ri = hardening.LinearIsotropicHardeningRule(sy_interp, hmodulus)
  hard_ri = hardening.Chaboche(iso_ri,
      [interpolate.ConstantInterpolate(0.0)], 
      [hardening.ConstantGamma(interpolate.ConstantInterpolate(0.0))],
      [interpolate.ConstantInterpolate(0.0)],
      [interpolate.ConstantInterpolate(1.0)])
  flow_ri = ri_flow.RateIndependentNonAssociativeHardening(surface, hard_ri)
  
  rate_independent = models.SmallStrainRateIndependentPlasticity(elastic_m,
      flow_ri)

  # Combined model
  combined = models.KMRegimeModel(elastic_m, [rate_independent, rate_dependent],
      [g0], kboltz, b, eps0)

  # Do things

  # Strain rate jump
  Ts = np.linspace(Tmin, Tmax, 5)
  erates = [1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6]
  es = [0.01] * len(erates)
  for T in Ts:
    strain, stress = rate_jump_test(combined, erates, es, T)
    plt.plot(strain, stress)
  
  legend_string = map(lambda x: "%3.0f K" % x, Ts)
  plt.legend(legend_string, loc = 'best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.title("Strain rate jump test")
  plt.show()

  # Temperature jump test (no thermal strains)
  Ts = list(np.linspace(Tmin, Tmax, 3))
  Ts += Ts[::-1]
  erates = [1.0e-6, 1.0e-4, 1.0e-2, 1.0]
  for er in erates:
    strain, stress = temp_jump_test(combined, er, es, Ts)
    plt.plot(strain, stress)
  
  legend_string = map(lambda x: "%1.0e /s" % x, erates)
  plt.legend(legend_string, loc = 'best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.title("Temperature jump test")
  plt.show()
