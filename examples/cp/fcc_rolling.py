#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal
from neml.math import rotations, tensors
from neml import elasticity

import matplotlib.pyplot as plt

if __name__ == "__main__":
  # Everyone's favorite FCC rolling example!
  N = 100
  
  nthreads = 4

  L = np.array([[0,0,0],[0,1.0,0],[0,0,-1.0]])
  erate = 1.0e-4
  steps = 100
  emax = 0.5

  E = 100000.0
  nu = 0.3

  t0 = 50.0
  ts = 50.0
  b = 100.0

  g0 = 1.0
  n = 12.0

  # Setup
  L *= erate
  dt = emax / steps / erate
  orientations = rotations.random_orientations(N)
  
  strengthmodel = slipharden.VoceSlipHardening(ts, b, t0)
  slipmodel = sliprules.PowerLawSlipRule(strengthmodel, g0, n)
  imodel = inelasticity.AsaroInelasticity(slipmodel)
  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  kmodel = kinematics.StandardKinematicModel(emodel, imodel)

  lattice = crystallography.CubicLattice(1.0)
  lattice.add_slip_system([1,1,0],[1,1,1])

  model = singlecrystal.SingleCrystalModel(kmodel, lattice)

  pmodel = polycrystal.TaylorModel(model, orientations)
  
  e = [0.0]
  s = [0.0]
  for i in range(steps):
    pmodel.deformation_step(L, dt, nthreads = nthreads)
    e.append(e[-1] + erate * dt)
    s.append(pmodel.s[2])
  
  plt.figure()
  plt.plot(e,s)
  plt.show()
