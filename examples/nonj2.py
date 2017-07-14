#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, neml, elasticity, drivers, surfaces, hardening, ri_flow

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  E = 150000.0
  nu = 0.3

  sY = 150.0
  h = 1.0e-2
  l = 1.0

  youngs = elasticity.YoungsModulus(E)
  poisson = elasticity.PoissonsRatio(nu)
  emodel = elasticity.IsotropicLinearElasticModel(youngs, poisson)

  surface = surfaces.IsoJ2I1(h, l)

  model = neml.SmallStrainPerfectPlasticity(emodel, surface, sY)

  res_tension = drivers.uniaxial_test(model, 1.0e-2)
  res_compres = drivers.uniaxial_test(model, 1.0e-2, 
      sdir = np.array([-1,0,0,0,0,0]))

  plt.plot(res_tension['strain'], res_tension['stress'], 'k-')
  plt.plot(res_compres['strain'], res_compres['stress'], 'r-')
  plt.show()
 

  E = 150000.0
  nu = 0.3

  sY = 150.0
  h = 1.0e-2
  l = 1.0
  K = E / 50.0

  youngs = elasticity.YoungsModulus(E)
  poisson = elasticity.PoissonsRatio(nu)
  emodel = elasticity.IsotropicLinearElasticModel(youngs, poisson)

  surface = surfaces.IsoKinJ2I1(h, l)
  
  hiso = hardening.LinearIsotropicHardeningRule(sY, -K/10)
  hkin = hardening.LinearKinematicHardeningRule(K)
  hardening = hardening.CombinedHardeningRule(hiso, hkin)
  flow = ri_flow.RateIndependentAssociativeFlow(surface, hardening)

  model = neml.SmallStrainRateIndependentPlasticity(emodel, flow)
  
  smax = 200.0
  Rs = [-1.0, -0.75, -0.5, -0.25, 0.0]
  srate = 1.0
  ncycles = 50
  for R in Rs:
    res = drivers.stress_cyclic(model, smax, R, srate, ncycles)
    plt.plot(res['cycles'], res['max'])

  plt.legend(map(str, Rs))

  plt.show()
