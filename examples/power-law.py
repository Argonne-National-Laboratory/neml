#!/usr/bin/env python3

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow, parse

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  E = 100000.0
  nu = 0.3
  s0 = 100.0
  A = 200.0
  n = 0.2

  elastic = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")

  surface = surfaces.IsoJ2()
  iso = hardening.PowerLawIsotropicHardeningRule(s0, A, n)

  flow = ri_flow.RateIndependentAssociativeFlow(surface, iso)
  model = models.SmallStrainRateIndependentPlasticity(elastic, flow)

  model2 = parse.parse_xml("example.xml", "powerlaw")
  
  erate = 1.0

  res = drivers.uniaxial_test(model, erate, emax = 0.1)
  res2 = drivers.uniaxial_test(model2, erate, emax = 0.1) 
  plt.plot(res['strain'], res['stress'])
  plt.plot(res2['strain'], res2['stress'], ls = '--')
  plt.show()
