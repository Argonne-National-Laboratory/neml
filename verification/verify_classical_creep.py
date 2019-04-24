#!/usr/bin/env python3

import sys
sys.path.append('..')

from neml import models, elasticity, drivers, damage, creep

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  # Elasticity
  E = 180000.0
  nu = 0.3

  # Power law creep parameters
  n = 6.0
  s0 = 100.0
  A = s0 ** (-n)
  A = 1.0e-22

  # Damage parameters
  xi = 0.478
  phi = 1.914
  S = 1.0e19

  # Loading
  srate = 100.0

  # Setup model
  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  bmodel = models.SmallStrainElasticity(emodel)
  scmodel = creep.PowerLawCreep(A, n)
  cfmodel = creep.J2CreepModel(scmodel)
  cmodel = models.SmallStrainCreepPlasticity(emodel, bmodel, cfmodel)
  model = damage.ClassicalCreepDamageModel_sd(emodel, S, xi, phi, cmodel)

  # Computed life
  srange = np.linspace(s0/2,s0*2, 10)
  tfs = S**(xi) / (1+phi) * srange**(-xi)
  
  slife = []
  for s,tf in zip(srange, tfs):
    res = drivers.creep(model, s, srate, tf * 1.25, verbose = False,
        logspace = False, check_dmg = True, nsteps = 1000)
    slife.append(res['rtime'][-1])
  
  plt.plot(srange, tfs/3600.0)
  plt.plot(srange, np.array(slife)/3600.0)
  plt.ylim([0,18000])
  plt.xlabel("Stress (MPa)")
  plt.ylabel("Failure time (hrs)")
  plt.show()
