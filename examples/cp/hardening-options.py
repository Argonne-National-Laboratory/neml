#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures
from neml.math import rotations, tensors, nemlmath, matrix
from neml import elasticity, drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 10
  nthreads = 2

  E = 160000.0
  nu = 0.31

  g0 = 1.0e-4
  n = 10.0

  K = E / 50.0
  s0 = 75.0 

  emax = 0.01
  erate = 1.0e-4
  R = -1
  ncycles = 3
  
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
  #     1) Initial strength of 37.5
  #     2) K = E / 50
  #     3) abs on the slip rates
  h2 = slipharden.GeneralLinearHardening(M, [s0/2] * lattice.ntotal, 
      absval = True)

  # Option 3: diagonal hardening with:
  #     1) Initial strength of 0
  #     2) K = E / 50
  #     3) abs on the slip rates
  h3 = slipharden.GeneralLinearHardening(M, [0] * lattice.ntotal, 
      absval = True) 

  # Option 4: diagonal hardening with:
  #     1) Initial strength of 0
  #     2) K = E / 50
  #     3) no abs on the slip rates
  h4 = slipharden.GeneralLinearHardening(M, [0] * lattice.ntotal, 
      absval = False)

  # Option 5: no hardening, strength of 75
  h5 = slipharden.FixedStrengthHardening([s0] * lattice.ntotal)

  # Option 6: no hardening, strength of 37.5
  h6 = slipharden.FixedStrengthHardening([s0/2] * lattice.ntotal)

  # Option 7: no hardening, strength of 0
  h7 = slipharden.FixedStrengthHardening([0] * lattice.ntotal)

  # The order is the backstrength, the isotropic strength, and the flow resistance
  strength775 = sliprules.KinematicPowerLawSlipRule(h7, h7, h5, g0, n)
  strength766 = sliprules.KinematicPowerLawSlipRule(h7, h6, h6, g0, n)

  strength771 = sliprules.KinematicPowerLawSlipRule(h7, h7, h1, g0, n)
  strength722 = sliprules.KinematicPowerLawSlipRule(h7, h2, h2, g0, n)
  strength735 = sliprules.KinematicPowerLawSlipRule(h7, h3, h5, g0, n)

  def make_model(smodel):
    imodel = inelasticity.AsaroInelasticity(smodel)
    emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
    kmodel = kinematics.StandardKinematicModel(emodel, imodel)
    model = singlecrystal.SingleCrystalModel(kmodel, lattice)
    return polycrystal.TaylorModel(model, orientations, nthreads = nthreads)
  
  """
  # Perfect plasticity options
  smodels = [strength775, strength766]
  names = ["775", "766"]

  models = [make_model(s) for s in smodels]

  results = [drivers.strain_cyclic(m, emax, R, erate, ncycles)
      for m in models]
  
  plt.figure()
  plt.title("Various perfect plasticity options")
  for res, name in zip(results, names):
    plt.plot(res['strain'], res['stress'], label = name)
  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.show()
  """

  # Pure isotropic options
  smodels = [strength771]
  names = ["771"]

  models = [make_model(s) for s in smodels]

  results = [drivers.strain_cyclic(m, emax, R, erate, ncycles, verbose = True)
      for m in models]
  
  plt.figure()
  plt.title("Various isotropic hardening options")
  for res, name in zip(results, names):
    plt.plot(res['strain'], res['stress'], label = name)
  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.show()
