#!/usr/bin/env python3

import sys
sys.path.append('../../..')

import numpy as np

from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal, crystaldamage
from neml.math import rotations, tensors, nemlmath
from neml import elasticity, drivers

import matplotlib.pyplot as plt

nthreads = 4



def load_unload(model, emax, erate, nsteps, E, ff = 1.0, verbose = False, miter = 50, T = 0):
  
  sdir = np.array([1.0,0,0,0,0,0])

  driver = drivers.Driver_sd(model, verbose = verbose, miter = miter)

  for e in np.linspace(0,emax,nsteps)[1:]:
    driver.erate_step(sdir, erate, e/erate, T)
  
  t = driver.t[-1]
  my_stress = driver.stress[-1][0]
  dt = my_stress / E / nsteps / erate * ff

  einc_guess = np.array([-1,0.5,0.5,0,0,0]) * 1.0e-4
  ainc_guess = -0.01

  for i in range(nsteps):
    t += dt
    if i == 0:
      driver.erate_step(-sdir, erate, t, T, einc_guess = einc_guess,
          ainc_guess = ainc_guess)
    else:
      driver.erate_step(-sdir, erate, t, T)

  return driver.strain[:,0], driver.stress[:,0]

if __name__ == "__main__":
  N = 50
  emax = 0.1
  erate = 1.0e-4
  nsteps = 100
  nthreads = 4

  orientations = rotations.random_orientations(N)

  lattice = crystallography.CubicLattice(1.0)
  lattice.add_slip_system([1,1,0],[1,1,1])

  nslip = lattice.ntotal

  t0 = 100.0
  b = 10.0
  sat = 50

  g0 = 1.0e-4
  n = 6.0

  mu = 29000.0
  E = 120000.0
  nu = 0.3

  c = 40
  beta = 3.0

  emodel = elasticity.CubicLinearElasticModel(E, nu, mu, "moduli")

  strengthmodel = slipharden.VoceSlipHardening(sat, b, t0)

  slipmodel = sliprules.PowerLawSlipRule(strengthmodel, g0, n)

  imodel = inelasticity.AsaroInelasticity(slipmodel)

  base_kin = kinematics.StandardKinematicModel(emodel, imodel)
  base_cp_model = singlecrystal.SingleCrystalModel(base_kin, lattice)
  base_model = polycrystal.TaylorModel(base_cp_model, orientations, 
      nthreads = nthreads)


  base_kin = kinematics.StandardKinematicModel(emodel, imodel)
  base_cp_model = singlecrystal.SingleCrystalModel(base_kin, lattice)
  base_model = polycrystal.TaylorModel(base_cp_model, orientations, 
      nthreads = nthreads)

  dmodel = crystaldamage.WorkPlaneDamage()
  func = crystaldamage.SigmoidTransformation(c, beta)
  dmg_fun = crystaldamage.PlanarDamageModel(dmodel, func, func, lattice)
  dmg_kin = kinematics.DamagedStandardKinematicModel(emodel, imodel, dmg_fun)
  dmg_cp_model = singlecrystal.SingleCrystalModel(dmg_kin, lattice)
  dmg_model = polycrystal.TaylorModel(dmg_cp_model, orientations,
      nthreads = nthreads)
  
  strain1, stress1 = load_unload(base_model, emax, erate, nsteps, E, 
      verbose = False, ff = 1.2)
  strain2, stress2 = load_unload(dmg_model, emax, erate, nsteps, E,
      verbose = False, ff = 1.4)
  
  plt.plot(strain1, stress1)
  plt.plot(strain2, stress2)
  plt.show()
