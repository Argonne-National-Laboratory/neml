#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow, visco_flow, general_flow

import matplotlib.pyplot as plt

import numpy as np

if __name__ == "__main__":
  E = 1000.0
  nu = 0.27

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  s0 = 10.0
  Kp = E / 10.0

  elastic = elasticity.IsotropicLinearElasticModel(mu, "shear", K, "bulk")
  model1 = models.SmallStrainElasticity(elastic)
  
  surface = surfaces.IsoJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)

  flow1 = ri_flow.RateIndependentAssociativeFlow(surface, iso)
  model2 = models.SmallStrainRateIndependentPlasticity(elastic, flow1)
  
  n = 10.0
  eta = 100.0
  g = visco_flow.GPowerLaw(n)
  vmodel1 = visco_flow.PerzynaFlowRule(surface, iso, g, eta)
  flow2 = general_flow.TVPFlowRule(elastic, vmodel1)
  model3 = models.GeneralIntegrator(elastic, flow2)
  
  mods = [model1, model2, model3]

  ey = s0 / E

  strain = ey
  stress = strain * E
  erate = 1.0e-2
  shoulde = stress * strain / 2.0
  shouldp = 0.0
  
  print("Elastic range:")
  for model in mods:
    res = drivers.uniaxial_test(model, erate, emax = strain, nsteps = 200)
    ef = res['energy_density'][-1]
    pf = res['plastic_work'][-1]
    print("\tStrain energy density: %f/%f" % (ef,shoulde))
    print("\tPlastic work: %f/%f" % (pf,shouldp))

  strain = 2 * ey
  erate = 1.0e-2
  
  dstrain = strain - ey

  shoulde = s0 * ey / 2.0 + dstrain * s0 + dstrain * dstrain * Kp / 2.0
  shouldp = dstrain * s0 + dstrain * dstrain * Kp / 2.0

  mods = [model2, model3]
  
  print("Plastic range:")
  for model in mods:
    res = drivers.uniaxial_test(model, erate, emax = strain)
    ef = res['energy_density'][-1]
    pf = res['plastic_work'][-1]
    print("\tStrain energy density: %f/%f" % (ef,shoulde))
    print("\tPlastic work: %f/%f" % (pf,shouldp))

  
  print("Plastic shakedown")
  iso2 = hardening.VoceIsotropicHardeningRule(s0, s0/5, s0)
  vmodel2 = visco_flow.PerzynaFlowRule(surface, iso2, g, eta)
  flow3 = general_flow.TVPFlowRule(elastic, vmodel2)
  model4 = models.GeneralIntegrator(elastic, flow3) 

  res = drivers.strain_cyclic(model4, 2*ey, -1.0, erate, 50)
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.show()

  plt.plot(res['cycles'][:-1], np.diff(res['energy_density']), 'k-')
  plt.plot(res['cycles'][:-1], np.diff(res['plastic_work']), 'r-')
  plt.show()
  
  print("Elastic cycles")
  res = drivers.stress_cyclic(model4, s0, -1.0, erate, 50)
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.show()
  
  plt.plot(res['cycles'], res['energy_density'], 'k-')
  plt.plot(res['cycles'], res['plastic_work'], 'r-')
  plt.show()

  print("Elastic shakedown")
  res = drivers.stress_cyclic(model3, 1.1*s0, -1.0, erate, 150)
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.show()
  
  plt.plot(res['cycles'], res['energy_density'], 'k-')
  plt.plot(res['cycles'], res['plastic_work'], 'r-')
  plt.show()
