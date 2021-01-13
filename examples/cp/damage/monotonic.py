#!/usr/bin/env python3

import sys
sys.path.append('../../..')

import numpy as np

from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal, crystaldamage
from neml.math import rotations, tensors, nemlmath
from neml import elasticity, drivers

import matplotlib.pyplot as plt

nthreads = 4

if __name__ == "__main__":
  N = 50
  emax = 0.122
  nsteps = 300
  erate = 1.0e-4
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

  res = drivers.uniaxial_test(base_model, erate, emax = emax, nsteps = nsteps,
      verbose = False)
  res2 = drivers.uniaxial_test(dmg_model, erate, emax = emax, nsteps = nsteps,
      verbose = True)
  
  plt.plot(res['strain'], res['stress'])
  plt.plot(res2['strain'], res2['stress'])
  plt.show()
