#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, neml, elasticity, drivers, surfaces, hardening, visco_flow, general_flow

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  E = 151000.0
  nu = 0.33

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  n = 12.0
  eta = 150.0
  k = 6.0
  C = 24800.0
  g0 = 300.0
  Q = 86 - k
  gs = 300.0
  b = 10.0
  beta = 0.0
  
  shear = elasticity.ShearModulus(mu)
  bulk = elasticity.BulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

  surface = surfaces.IsoKinJ2()
  iso = hardening.VoceIsotropicHardeningRule(k, Q, b)
  cs = [C]
  gs = [hardening.SatGamma(gs, g0, beta)]
  hmodel = hardening.Chaboche(iso, cs, gs)

  fluidity = visco_flow.ConstantFluidity(eta)

  vmodel = visco_flow.ChabocheFlowRule(
      surface, hmodel, fluidity, n)

  flow = general_flow.TVPFlowRule(elastic, vmodel)

  model = neml.GeneralIntegrator(flow, verbose = False)
  
  # Uniaxial stress/strain curves at decades of strain rates
  erates = np.logspace(-6,2,9)
  for rate in erates:
    res = drivers.uniaxial_test(model, rate, verbose = False)
    plt.plot(res['strain'], res['stress'])
  
  plt.xlabel("Strain (-/-)")
  plt.ylabel("Stress (MPa)")
  plt.show()
  
  # A strain-controlled cyclic test
  res = drivers.strain_cyclic(model, 0.001, -0.25, 1.0e-4, 50,
      verbose = False)
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.xlabel("Strain (-/-)")
  plt.ylabel("Stress (MPa)")
  plt.show()

  plt.plot(res['cycles'], res['max'], 'k-')
  plt.plot(res['cycles'], res['min'], 'k-')
  plt.plot(res['cycles'], res['mean'], 'k-')
  plt.xlabel("Cycle")
  plt.ylabel("Stress (MPa)")
  plt.show()
  
  # A stress-controlled cyclic test
  res = drivers.stress_cyclic(model, 100.0, -1.0, 5.0, 50,
      hold_time = 10, verbose = False)
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.xlabel("Strain (-/-)")
  plt.ylabel("Stress (MPa)")
  plt.show()

  plt.plot(res['cycles'], res['max'], 'k-')
  plt.plot(res['cycles'], res['min'], 'k-')
  plt.plot(res['cycles'], res['mean'], 'k-')
  plt.xlabel("Cycle")
  plt.ylabel("Strain (-/-)")
  plt.show()
  

  # Stress relaxation test
  res = drivers.stress_relaxation(model, 0.02, 1.0e-4, 2000.0, verbose = False)
  plt.plot(res['time'], res['stress'], 'k-')
  plt.xlabel('Time (s)')
  plt.ylabel('Stress (MPa)')
  plt.show()

  plt.plot(res['rtime'], res['rrate'], 'k-')
  plt.xlabel('Time (s)')
  plt.ylabel('Relaxation rate (MPa/s)')
  plt.show()

  # Creep test
  res = drivers.creep(model, 150.0, 10.0, 6000.0)
  plt.plot(res['time'], res['strain'], 'k-')
  plt.xlabel('Time (s)')
  plt.ylabel('Strain (-/-)')
  plt.show()

  plt.plot(res['rtime'], res['rrate'], 'k-')
  plt.xlabel('Time (s)')
  plt.ylabel('Creep rate (1/s)')
  plt.show()
