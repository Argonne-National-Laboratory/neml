#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  E = 200000.0
  nu = 0.27

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  s0 = 300.0
  Kp = 0.0
  c = [30000.0]
  r = [60.0]
  
  elastic = elasticity.IsotropicLinearElasticModel(mu, "shear", K, "bulk")

  surface = surfaces.IsoKinJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
  gmodels = [hardening.ConstantGamma(g) for g in r]
  As = [0.0]
  ns = [1.0]
  hrule = hardening.Chaboche(iso, c, gmodels, As, ns)

  flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hrule)
  model = models.SmallStrainRateIndependentPlasticity(elastic, flow, verbose = False,
      check_kt = False)
  
  
  # Uniaxial stress/strain curves at decades of strain rates
  erates = np.logspace(-6,2,9)
  for rate in erates:
    res = drivers.uniaxial_test(model, rate)
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
