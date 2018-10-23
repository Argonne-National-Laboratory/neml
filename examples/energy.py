#!/usr/bin/env python

import sys
sys.path.append('..')

import numpy as np
import scipy.integrate as quad

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow, creep, visco_flow, general_flow

import matplotlib.pyplot as plt

if __name__ == "__main__":
  E = 150000.0
  nu = 0.3
  sY = 200.0
  H = 0.0001

  elastic = elasticity.IsotropicLinearElasticModel(E, "youngs",
      nu, "poissons")
  surface = surfaces.IsoJ2()
  iso = hardening.LinearIsotropicHardeningRule(sY, H)
  flow = ri_flow.RateIndependentAssociativeFlow(surface, iso)

  model1 = models.SmallStrainRateIndependentPlasticity(elastic, flow)

  model2 = models.SmallStrainPerfectPlasticity(elastic, surface, sY)
 
  n = 45.0
  eta = 200.0
  iso2 = hardening.LinearIsotropicHardeningRule(0, H)
  gpower = visco_flow.GPowerLaw(n, eta)
  vflow = visco_flow.PerzynaFlowRule(surface, iso2, gpower)
  gflow = general_flow.TVPFlowRule(elastic, vflow)
  model3 = models.GeneralIntegrator(elastic, gflow)

  # T is in hours, strain in percent, stress in MPa
  A = 1.85e-9
  n = 2.5
  m = 0.3

  smodel = creep.NortonBaileyCreep(A, n, m)
  cmodel = creep.J2CreepModel(smodel)
  model4 = models.SmallStrainCreepPlasticity(elastic, model1, cmodel)

  erate = 1.0e-4
  emax = 0.025

  energy = emax * sY - 0.5 * sY * sY/E

  dissipated = energy - 0.5 * sY * sY / E
  
  res1 = drivers.uniaxial_test(model1, erate, emax = emax)
  res2 = drivers.uniaxial_test(model2, erate, emax = emax)
  res3 = drivers.uniaxial_test(model3, erate, emax = emax)
  res4 = drivers.uniaxial_test(model4, erate, emax = emax)
  
  print("")

  print("Total energy")
  print("------------")
  print("Analytic:\t\t%f" % energy)
  print("Rate independent:\t%f" % res1['energy_density'][-1])
  print("Perfectly plastic:\t%f" % res2['energy_density'][-1])
  print("Rate dependent:\t\t%f*" % res3['energy_density'][-1])
  print("Creep + plasticity:\t%f*" % res4['energy_density'][-1])
  print("")

  print("Dissipation")
  print("------------")
  print("Analytic:\t\t%f" % dissipated)
  print("Rate independent:\t%f" % res1['plastic_work'][-1])
  print("Perfectly plastic:\t%f" % res2['plastic_work'][-1])
  print("Rate dependent:\t\t%f*" % res3['plastic_work'][-1])
  print("Creep + plasticity:\t%f*" % res4['plastic_work'][-1])
  print("")
  print("*: should not be exact")

  plt.plot(res1['strain'], res1['stress'], 'k-')
  plt.plot(res2['strain'], res2['stress'], 'k--')
  plt.plot(res3['strain'], res3['stress'], 'k:')
  plt.plot(res4['strain'], res4['stress'], 'k-.')
  plt.show()

  plt.plot(res1['strain'], res1['energy_density'], 'k-')
  plt.plot(res1['strain'], res1['plastic_work'], 'r-')

  plt.plot(res2['strain'], res2['energy_density'], 'k--')
  plt.plot(res2['strain'], res2['plastic_work'], 'r--')

  plt.plot(res3['strain'], res3['energy_density'], 'k:')
  plt.plot(res3['strain'], res3['plastic_work'], 'r:')

  plt.plot(res4['strain'], res4['energy_density'], 'k-.')
  plt.plot(res4['strain'], res4['plastic_work'], 'r-.')

  plt.show()
