#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal
from neml.math import rotations, tensors, nemlmath
from neml import elasticity

import matplotlib.pyplot as plt

if __name__ == "__main__":
  sdir = np.array([0,0,1.0,0,0,0])
  erate = 1.0e-4
  steps = 200
  emax = 0.5

  E = 146000.0
  nu = 0.3

  mu =  E / (2 * (1 + nu))

  alpha = 0.3
  burg = 20.2e-9
  
  rhos = [9e10,1e11,2e11,3e11]
  
  strains = []
  stresses = []

  for rho in rhos:
    t0 = alpha*mu*burg*np.sqrt(rho)
    print("HMM",t0)
    ts = 100.0
    b = 0.4

    g0 = 1.0
    n = 6.7

    N = 100

    nthreads = 20

    strengthmodel = slipharden.VoceSlipHardening(ts, b, t0)
    slipmodel = sliprules.PowerLawSlipRule(strengthmodel, g0, n)
    imodel = inelasticity.AsaroInelasticity(slipmodel)
    emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
    kmodel = kinematics.StandardKinematicModel(emodel, imodel)

    lattice = crystallography.CubicLattice(1.0)
    lattice.add_slip_system([1,1,0],[1,1,1])

    model = singlecrystal.SingleCrystalModel(kmodel, lattice, verbose = False)

    orientations = rotations.random_orientations(N)

    dt = emax / erate / steps

    pmodel = polycrystal.TaylorModel(model, orientations)

    e = [0.0]
    s = [0.0]
    for i in range(steps):
      print(i)
      pmodel.take_stress_step(sdir, erate, dt, nthreads = nthreads, 
          verbose = True)
      e.append(pmodel.e[2])
      s.append(pmodel.s[2])

    strains.append(e)
    stresses.append(s)
  
  plt.style.use("presentation")
  plt.figure()
  for e,s,r in zip(strains, stresses, rhos):
    plt.plot(e,s, label = r"$\rho=%3.1e\,\mathrm{m^{-2}}$" % r)

  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.legend(loc='best')
  plt.tight_layout()
  plt.show()
