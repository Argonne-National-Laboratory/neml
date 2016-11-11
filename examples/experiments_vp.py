#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import neml, elasticity, drivers, surfaces, hardening, visco_flow

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  E = 200000.0
  nu = 0.27

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  s0 = 86.0
  Kp = E/100.0

  n = 2.0
  eta = 150.0
  
  shear = elasticity.ConstantShearModulus(mu)
  bulk = elasticity.ConstantBulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)
  surface = surfaces.IsoJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
  g = visco_flow.GPowerLaw(n)

  flow = visco_flow.PerzynaFlowRule(surface, iso, g, eta)
  model = neml.SmallStrainViscoPlasticity(elastic, flow, verbose = True)
  
  # Uniaxial stress/strain curves at decades of strain rates
  erates = np.logspace(-6,2,9)
  for rate in erates:
    res = drivers.uniaxial_test(model, rate, verbose = True)
    plt.plot(res['strain'], res['stress'])
  
  plt.xlabel("Strain (-/-)")
  plt.ylabel("Stress (MPa)")
  plt.show()
  
  # A strain-controlled cyclic test
  res = drivers.strain_cyclic(model, 0.01, -0.25, 1.0e-4, 50)
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
  res = drivers.stress_cyclic(model, 525.0, -1.0, 5.0, 50,
      hold_time = 100)
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
  res = drivers.stress_relaxation(model, 0.02, 1.0e-4, 2000.0)
  plt.plot(res['time'], res['stress'], 'k-')
  plt.xlabel('Time (s)')
  plt.ylabel('Stress (MPa)')
  plt.show()

  plt.plot(res['rtime'], res['rrate'], 'k-')
  plt.xlabel('Time (s)')
  plt.ylabel('Relaxation rate (MPa/s)')
  plt.show()

  # Creep test
  res = drivers.creep(model, 350.0, 10.0, 2000.0)
  plt.plot(res['time'], res['strain'], 'k-')
  plt.xlabel('Time (s)')
  plt.ylabel('Strain (-/-)')
  plt.show()

  plt.plot(res['rtime'], res['rrate'], 'k-')
  plt.xlabel('Time (s)')
  plt.ylabel('Creep rate (1/s)')
  plt.show()
