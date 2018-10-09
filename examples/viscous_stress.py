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
  R = 100.0
  d = 10.0
  n = 3.0
  eta = 5000.0

  elastic = elasticity.IsotropicLinearElasticModel(E, "youngs",
      nu, "poissons")
  surface = surfaces.IsoJ2()
  iso = hardening.VoceIsotropicHardeningRule(sY, R, d)
  flow = ri_flow.RateIndependentAssociativeFlow(surface, iso)

  ri_model = models.SmallStrainRateIndependentPlasticity(elastic, flow)

  gpower = visco_flow.GPowerLaw(n)
  vflow = visco_flow.PerzynaFlowRule(surface, iso, gpower, eta)
  gflow = general_flow.TVPFlowRule(elastic, vflow)
  rd_model = models.GeneralIntegrator(elastic, gflow)


  eps_dot = 1.0e2
  emax = 0.5

  res_ri = drivers.uniaxial_test(ri_model, eps_dot, emax = 0.5)
  res_rd = drivers.uniaxial_test(rd_model, eps_dot, emax = 0.5)

  plt.plot(res_ri['strain'], res_ri['stress'], 'k-')
  plt.plot(res_rd['strain'], res_rd['stress'], 'r-')
  
  plt.plot(res_ri['strain'], res_rd['stress'] - res_ri['stress'], 'b--')
  
  sigma_v = np.sqrt(3.0 / 2.0) * (eta*eps_dot) ** (1.0 / n)

  plt.plot(res_ri['strain'], np.ones(res_ri['strain'].shape) * sigma_v, 'b:')

  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")

  plt.legend(['Rate independent', 'Rate dependent', 'Difference', 
    'Viscous stress'], loc = 'best')

  plt.show()

