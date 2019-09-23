#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures
from neml.math import rotations, tensors
from neml import elasticity

import matplotlib.pyplot as plt

if __name__ == "__main__":
  # Everyone's favorite FCC rolling example!
  N = 200
  
  nthreads = 4

  L = np.array([[0,0,0],[0,1.0,0],[0,0,-1.0]])
  erate = 1.0e-4
  steps = 100
  emax = 1.0

  E = 100000.0
  nu = 0.3

  t0 = 30.0
  ts = 10.0
  b = 1.0

  g0 = 1.0
  n = 12.0

  # Setup
  L *= erate
  dt = emax / steps / erate
  orientations = rotations.random_orientations(N)
 
  lattice = crystallography.CubicLattice(1.0)
  lattice.add_slip_system([1,1,0],[1,1,1])

  strengthmodel = slipharden.VoceSlipHardening(ts, b, t0)
  slipmodel = sliprules.PowerLawSlipRule(strengthmodel, g0, n)
  imodel = inelasticity.AsaroInelasticity(slipmodel)
  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  kmodel = kinematics.StandardKinematicModel(emodel, imodel)

  model = singlecrystal.SingleCrystalModel(kmodel, lattice)

  pmodel = polycrystal.TaylorModel(model, orientations)
 
  e = [0.0]
  s = [0.0]
  for i in range(steps):
    print(i)
    pmodel.deformation_step(L, dt, nthreads = nthreads)
    e.append(e[-1] + erate * dt)
    s.append(pmodel.s[2])
  
  polefigures.pole_figure_discrete(orientations,[1,1,1],lattice)
  plt.title("Initial, <111>")
  plt.show()

  polefigures.pole_figure_discrete(pmodel.orientations,[1,1,1],lattice)
  plt.title("Final, <111>")
  plt.show()

  polefigures.pole_figure_discrete(orientations,[1,0,0],lattice)
  plt.title("Initial, <100>")
  plt.show()

  polefigures.pole_figure_discrete(pmodel.orientations,[1,0,0],lattice)
  plt.title("Final, <100>")
  plt.show()
