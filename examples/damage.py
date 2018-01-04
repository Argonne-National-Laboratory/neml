#!/usr/bin/env python

import sys
sys.path.append('..')

import numpy as np
from neml import interpolate, solvers, neml, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep, damage, drivers

import matplotlib.pyplot as plt

def simple_ex():
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

def creep_ex():
  E = 92000.0
  nu = 0.3

  s0 = 120.0

  A = 1.0e-10
  n = 3.0

  Kp = E/500
  H = E/500

  smodel = creep.PowerLawCreep(A, n)
  cmodel = creep.J2CreepModel(smodel)

  E = elasticity.YoungsModulus(E)
  v = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(E, v)

  surface = surfaces.IsoKinJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
  kin = hardening.LinearKinematicHardeningRule(H)
  hrule = hardening.CombinedHardeningRule(iso, kin)

  flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

  pmodel = neml.SmallStrainRateIndependentPlasticity(elastic, 
      flow)
  
  bmodel = neml.SmallStrainCreepPlasticity(pmodel, cmodel)

  A_damg = 1.0e-2
  a_damg = 1.0

  model = damage.NEMLPowerLawDamagedModel_sd(A_damg, a_damg, bmodel, elastic,
      verbose = False)

  #res = drivers.uniaxial_test(model, 1.0e-2, emax = 0.25)
  #plt.plot(res['strain'], res['stress'])
  #plt.show()

  res = drivers.creep(model, 120.0, 1.0, 393.0, verbose = False)
  
  plt.plot(res['rtime'], res['rstrain'])
  plt.show()
  plt.loglog(res['rtime'], res['rrate'])
  plt.show()


if __name__ == "__main__":
  #simple()
  creep_ex()
