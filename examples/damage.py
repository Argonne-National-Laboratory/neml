#!/usr/bin/env python

import sys
sys.path.append('..')

import numpy as np
from neml import interpolate, solvers, neml, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep, damage, drivers

import matplotlib.pyplot as plt

def simple():
  E = 92000.0
  nu = 0.3

  s0 = 180.0
  Kp = 1000.0
  H = 1000.0

  E = elasticity.YoungsModulus(E)
  v = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(E, v)

  surface = surfaces.IsoKinJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
  kin = hardening.LinearKinematicHardeningRule(H)
  hrule = hardening.CombinedHardeningRule(iso, kin)

  flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

  bmodel = neml.SmallStrainRateIndependentPlasticity(elastic, 
      flow)

  A = 0.0e-6
  a = 2.2
  model_off = damage.NEMLPowerLawDamagedModel_sd(A, a, bmodel, elastic)

  A = 2e-5
  model_on = damage.NEMLPowerLawDamagedModel_sd(A, a, bmodel, elastic)

  res_off = drivers.uniaxial_test(model_off, 1.0e-2, emax = 0.13)
  res_on = drivers.uniaxial_test(model_on, 1.0e-2, emax = 0.13)

  plt.plot(res_off['strain'], res_off['stress'], 'k-')
  plt.plot(res_on['strain'], res_on['stress'], 'r-')
  plt.show()

if __name__ == "__main__":
  simple()
