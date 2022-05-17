#!/usr/bin/env python3

import sys
sys.path.append('..')

from neml import models, elasticity, drivers, damage, creep

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  s = 150.0

  # We'll do it in terms of Hayhurst-Leckie
  nu = 1.8
  eta = 2.1
  s0 = 7000.0
  w0 = 3.25e-2
  n = 4.0
  e0 = 1.0 # This is key

  # My parameters
  xi = nu
  phi = eta
  A = s0 / (w0**(1.0/nu))

  tf = 10000

  # Hayhurst solution
  times = np.linspace(0,tf,100)
  dmg = (1-(eta+1)*w0*(s/s0)**nu*times)**(1.0/(eta+1.0))
  strain = -e0 * (s/s0)**(-nu) * dmg**(eta+1) / ((1 + eta - n) * w0) * (s / (s0*dmg)
      )**n + e0 * (s/s0)**(n-nu)/((1+eta-n)*w0)

  
  E = 10000000.0
  nu = 0.3

  srate = 100.0

  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  bmodel = models.SmallStrainElasticity(emodel)
  scmodel = creep.NormalizedPowerLawCreep(s0, n)
  cfmodel = creep.J2CreepModel(scmodel)
  cmodel = models.SmallStrainCreepPlasticity(emodel, bmodel, cfmodel)
  model = damage.NEMLScalarDamagedModel_sd(emodel, cmodel,
          damage.ModularCreepDamage(emodel, A, xi, phi,
              damage.VonMisesEffectiveStress()))

  res = drivers.creep(model, s, srate, tf, nsteps = 1000)

  plt.plot(times, strain)
  plt.plot(res['rtime'], res['rstrain'])
  plt.show()
