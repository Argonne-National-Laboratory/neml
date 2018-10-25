#!/usr/bin/env python

import sys
sys.path.append('..')

import numpy as np
from neml import interpolate, solvers, models, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep, damage, drivers

import matplotlib.pyplot as plt

def simple_ex():
  E = 92000.0
  nu = 0.3

  s0 = 180.0
  Kp = 1000.0
  H = 1000.0

  elastic = elasticity.IsotropicLinearElasticModel(E, "youngs",
      nu, "poissons")

  surface = surfaces.IsoKinJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
  kin = hardening.LinearKinematicHardeningRule(H)
  hrule = hardening.CombinedHardeningRule(iso, kin)

  flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

  bmodel = models.SmallStrainRateIndependentPlasticity(elastic, 
      flow)

  A = 0.0e-6
  a = 2.2
  model_off = damage.NEMLPowerLawDamagedModel_sd(elastic, A, a, bmodel)

  A = 2e-5
  model_on = damage.NEMLPowerLawDamagedModel_sd(elastic, A, a, bmodel)

  res_off = drivers.uniaxial_test(model_off, 1.0e-2, emax = 0.13)
  res_on = drivers.uniaxial_test(model_on, 1.0e-2, emax = 0.13)

  plt.plot(res_off['strain'], res_off['stress'], 'k-')
  plt.plot(res_on['strain'], res_on['stress'], 'r-')
  plt.show()

def unload_ex():
  E = 92000.0
  nu = 0.3

  s0 = 180.0
  Kp = 1000.0
  H = 1000.0

  elastic = elasticity.IsotropicLinearElasticModel(E, "youngs",
      nu, "poissons")

  surface = surfaces.IsoKinJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
  kin = hardening.LinearKinematicHardeningRule(H)
  hrule = hardening.CombinedHardeningRule(iso, kin)

  flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

  bmodel = models.SmallStrainRateIndependentPlasticity(elastic, 
      flow)

  A = 2e-5
  a = 2.2
  model = damage.NEMLPowerLawDamagedModel_sd(elastic, A, a, bmodel)

  driver = drivers.Driver_sd(model)
  nsteps = 25
  sdir = np.array([1,0,0,0,0,0])
  erate = 1.0e-5
  e_inc = 0.1 / nsteps
  for i in range(nsteps):
    einc, ainc = driver.erate_einc_step(sdir, erate, e_inc, 0.0)

  estrain = model.elastic_strains(driver.stress_int[-1], driver.T_int[-1],
      driver.stored[-1])

  print("Calculated elastic strain: %f" % estrain[0])

  nsteps = 20
  dt = 0.1
  rstress = driver.stress_int[-1]
  rstrain = driver.strain_int[-1][0]
  for m in np.linspace(0,1.0,nsteps, endpoint = False)[::-1]:
    driver.stress_step(rstress*m, driver.t_int[-1] + dt, driver.T_int[-1])

  fstrain = driver.strain_int[-1][0]

  print("Actual elastic strain: %f" % (rstrain-fstrain))

  print("Calculated final elastic strain: %f" % model.elastic_strains(
    driver.stress_int[-1], driver.T_int[-1], driver.stored[-1])[0])

  plt.plot(driver.strain[:,0], driver.stress[:,0])
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

  elastic = elasticity.IsotropicLinearElasticModel(E, "youngs",
      nu, "poissons")

  surface = surfaces.IsoKinJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
  kin = hardening.LinearKinematicHardeningRule(H)
  hrule = hardening.CombinedHardeningRule(iso, kin)

  flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)

  pmodel = models.SmallStrainRateIndependentPlasticity(elastic, 
      flow)
  
  bmodel = models.SmallStrainCreepPlasticity(elastic, pmodel, cmodel)

  A_damg = 1.0e-2
  a_damg = 1.0

  model = damage.NEMLPowerLawDamagedModel_sd(elastic, A_damg, a_damg, bmodel,
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
  simple_ex()
  creep_ex()
  unload_ex()
