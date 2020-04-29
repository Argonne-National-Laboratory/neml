#!/usr/bin/env python3

import sys
sys.path.append('..')

import numpy as np
from neml import interpolate, solvers, models, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep, damage, drivers

import matplotlib.pyplot as plt

def workrate_ex():
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

  A_d   = 4.0
  P_d   = 1e4
  n_d   = 2.0
  workrate = interpolate.WorkRateFunc(A_d,P_d,n_d)

  Q = 1e3
  m = 2.0
  G = 1e-9
  H = 1.0
  xi = 1.0
  phi = 1.0
  model = damage.NEMLWorkRateFunctionDamage_sd(elastic, workrate, Q, m, G, H,
                                              xi, phi, bmodel)

  res = drivers.uniaxial_test(model, 1.0e-3, emax = 0.03)
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.show()

  res = drivers.creep(model, 100.0, 5.0, 1e6, check_dmg = True)
  plt.semilogx(res['rtime'], res['rstrain'], 'k-')
  plt.show()


if __name__ == "__main__":
  workrate_ex()
