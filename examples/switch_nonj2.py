#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, neml, elasticity, drivers, surfaces, hardening, ri_flow, interpolate, visco_flow, general_flow

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  # Model meta parameters
  Tmin = 290.0
  Tmax = 873.0
  nrate_sens = 40

  # Fixed parameters
  E_poly = [-78.275924, 236950.5588]
  nu = 0.3

  b = 0.248 * 1.0e-6
  kboltz = 1.38064e-23 * 1000.0
  eps0 = 1.0e10

  A = -9.774577
  B = -1.407274
  C = -5.060187
  g0 = 0.373716

  K = 150000.0 / 50.0
  h = 1.0e-2
  l = 1.0

  # Elastic
  E_m = elasticity.YoungsModulus(interpolate.PolynomialInterpolate(E_poly))
  nu_m = elasticity.PoissonsRatio(nu)
  elastic_m = elasticity.IsotropicLinearElasticModel(E_m, nu_m)

  # Rate sensitivity interpolates
  Ts_rate = np.linspace(Tmin, Tmax, num = nrate_sens)
  n_values = []
  eta_values = []
  flow_stress_values = []
  for T in Ts_rate:
    mu = elastic_m.G(T)
    n_values.append(-mu*b**3.0 / (kboltz * T * A))
    eta_values.append(np.exp(B) * eps0 ** (kboltz * T * A / (mu * b**3.0)) * mu)
    flow_stress_values.append(mu * np.exp(C))

  n_interp = interpolate.PiecewiseLinearInterpolate(list(Ts_rate), n_values)
  eta_interp = interpolate.PiecewiseLinearInterpolate(list(Ts_rate), eta_values)
  eta_m = visco_flow.ConstantFluidity(eta_interp)
  flow_interp = interpolate.PiecewiseLinearInterpolate(list(Ts_rate), flow_stress_values)
  
  iso_rd = hardening.LinearIsotropicHardeningRule(0.0, -K/10)
  iso_ri = hardening.LinearIsotropicHardeningRule(flow_interp, 
      interpolate.ConstantInterpolate(-K/10))
  
  hmodel_rd = hardening.Chaboche(iso_rd, [K * 3.0/2.0], [hardening.ConstantGamma(0.0)])
  hmodel_ri = hardening.Chaboche(iso_ri, [K * 3.0/2.0], [hardening.ConstantGamma(0.0)])

  surface_m = surfaces.IsoKinJ2()

  visco_flow_m = visco_flow.ChabocheFlowRule(surface_m, hmodel_rd, eta_m, 
      n_interp)
  rd_flow = general_flow.TVPFlowRule(elastic_m, visco_flow_m)
  rd_model = neml.GeneralIntegrator(elastic_m, rd_flow)
  
  ri_flow_m = ri_flow.RateIndependentNonAssociativeHardening(surface_m, 
      hmodel_ri)
  ri_model = neml.SmallStrainRateIndependentPlasticity(elastic_m,
      ri_flow_m)

  model = neml.KMRegimeModel(elastic_m, [ri_model, rd_model], [g0],
      kboltz, b, eps0)

  smax = 400.0
  Rs = [-1.0, -0.75, -0.5, -0.25, 0.0]
  srate = 1.0e-3
  ncycles = 50
  for R in Rs:
    res = drivers.stress_cyclic(model, smax, R, srate, ncycles, T = 550+273.15)
    plt.plot(res['cycles'], res['max'])
    #plt.plot(res['strain'], res['stress'])
    #plt.show()


  plt.legend(map(str, Rs))

  plt.show()
