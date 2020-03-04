#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures
from neml.math import rotations, tensors, nemlmath, matrix
from neml import elasticity, drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 100
  nthreads = 1

  E = 160000.0
  nu = 0.31

  g0 = 1.0
  n = 10.0

  K = E / 50.0
  s0 = 75.0 

  emax = 0.01
  erate = 1.0e-4
  R = -1
  ncycles = 2
  
  lattice = crystallography.CubicLattice(1.0)
  lattice.add_slip_system([1,1,0],[1,1,1])

  orientations = rotations.random_orientations(N)

  # Matrix used in the hardening options
  M = matrix.SquareMatrix(lattice.ntotal, type = "diagonal", 
      data = [K] * lattice.ntotal)

  # Option 1: diagonal hardening with:
  #     1) Initial strength of 75
  #     2) K = E / 50
  #     3) abs on the slip rates
  h1 = slipharden.GeneralLinearHardening(M, [s0] * lattice.ntotal, 
      absval = True)

  # Option 2: diagonal hardening with:
  #     1) Initial strength of 0
  #     2) K = E / 50
  #     3) no abs on the slip rates
  h2 = slipharden.GeneralLinearHardening(M, [0] * lattice.ntotal, 
      absval = False)

  # Option 3: no hardening, strength of 75
  h3 = slipharden.FixedStrengthHardening([s0] * lattice.ntotal)

  # Option 4: no hardening, strength of 0
  h4 = slipharden.FixedStrengthHardening([0] * lattice.ntotal)

  # The order is the backstrength, the isotropic strength, and the flow resistance
  strength441 = sliprules.PowerLawSlipRule(h3, g0, n)

  # Get the actual polycrystal models
  smodels = [strength441]
  names = ["441"]

  def make_model(smodel):
    imodel = inelasticity.AsaroInelasticity(smodel)
    emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
    kmodel = kinematics.StandardKinematicModel(emodel, imodel)
    model = singlecrystal.SingleCrystalModel(kmodel, lattice, verbose = True)
    return model
    return polycrystal.TaylorModel(model, orientations, nthreads = nthreads)

  models = [make_model(s) for s in smodels]

  results = [drivers.uniaxial_test(m, erate, verbose = True) for m in models]

  plt.plot(results[0]['strain'], results[0]['stress'])
  plt.show()
