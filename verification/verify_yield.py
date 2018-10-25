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
  Kp = E / 100.0

  elastic = elasticity.IsotropicLinearElasticModel(mu, "shear", K, "bulk")
  surface = surfaces.IsoJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)

  flow = ri_flow.RateIndependentAssociativeFlow(surface, iso)

  model = models.SmallStrainRateIndependentPlasticity(elastic, flow)

  erate = 1.0e-4
  res = drivers.uniaxial_test(model, erate, nsteps = 1000)
  strain = res['strain']
  stress = res['stress']


  print("Young's modulus: %f / %f" % (E, stress[1]/strain[1]))
  
  ys_tol = 1.0e-6
  tols = np.abs(stress - strain*E)
  ysi = np.where(tols>ys_tol)[0][0]
  ys = stress[ysi]
  
  print("Yield stress: %f / %f" % (s0, ys))

  print("Hardening modulus: %f / %f" % (Kp, (stress[ysi+2]-stress[ysi+1])/(strain[ysi+2]-strain[ysi+1])))

  plt.plot(res['strain'], res['stress'])
  plt.show()

