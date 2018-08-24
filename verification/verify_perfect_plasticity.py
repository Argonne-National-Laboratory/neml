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

  shear = elasticity.ShearModulus(mu)
  bulk = elasticity.BulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

  surface = surfaces.IsoJ2()

  model = models.SmallStrainPerfectPlasticity(elastic, surface, s0)
  
  erate = 1.0e-2
  strain = 0.02
  res = drivers.uniaxial_test(model, erate, emax = strain)

  ye = s0 / E
  energy = 0.5*ye*s0 + (strain - ye) * s0
  work = (strain - ye) * s0

  print("Strain energy: %f/%f" % (res['energy_density'][-1], energy))
  print("Plastic work: %f/%f" % (res['plastic_work'][-1], work))

  plt.plot(res['strain'], res['stress'], 'k-')
  plt.show()

  plt.plot(res['strain'], res['energy_density'], 'k-')
  plt.plot(res['strain'], res['plastic_work'], 'r-')
  plt.show()
