#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow, interpolate

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  E = 200000.0
  nu = 0.27

  points = [0.0, 0.01, 0.02, 0.045]
  values = [100.0, 150.0, 160.0, 170.0]
  youngs = elasticity.YoungsModulus(E)
  poissons = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(youngs, poissons)
  surface = surfaces.IsoJ2()
  ifn = interpolate.PiecewiseLinearInterpolate(points, values)
  iso = hardening.InterpolatedIsotropicHardeningRule(ifn)
  flow = ri_flow.RateIndependentAssociativeFlow(surface, iso)
  model = models.SmallStrainRateIndependentPlasticity(elastic, flow)

  res = drivers.uniaxial_test(model, 1.0e-4, emax = 0.05)
  
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.plot(res['strain'] + values[0] / E, [ifn(e) for e in res['strain']], 'r-')
  plt.xlim([0,0.05])
  plt.ylim([0,175])
  plt.show()
