#!/usr/bin/env python3

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  E = 150000.0
  nu = 0.3

  sY = 150.0
  h = 1.0e-2
  l = 1.0

  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")

  surface = surfaces.IsoJ2I1(h, l)

  model = models.SmallStrainPerfectPlasticity(emodel, surface, sY)

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

  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")

  surface = surfaces.IsoKinJ2I1(h, l)
  
  hiso = hardening.LinearIsotropicHardeningRule(sY, -K/10)
  
  hmodel = hardening.Chaboche(hiso, [K * 3.0/2.0], 
      [hardening.ConstantGamma(0.0)], [0.0], [1.0])
  flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hmodel)

  #hkin = hardening.LinearKinematicHardeningRule(K)
  #hmodel = hardening.CombinedHardeningRule(hiso, hkin)
  #flow = ri_flow.RateIndependentAssociativeFlow(surface, hmodel)

  model = models.SmallStrainRateIndependentPlasticity(emodel, flow)
  
  smax = 200.0
  Rs = [-1.0, -0.75, -0.5, -0.25, 0.0]
  srate = 1.0
  ncycles = 50
  for R in Rs:
    res = drivers.stress_cyclic(model, smax, R, srate, ncycles)
    plt.plot(res['cycles'], res['max'])

  plt.legend(map(str, Rs))

  plt.show()
