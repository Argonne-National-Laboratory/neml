#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures
from neml.math import rotations, tensors, nemlmath
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

  T = 0

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

  pmodel = polycrystal.TaylorModel(model, orientations, nthreads = nthreads)

  h_n = pmodel.init_store()
 
  d_inc = nemlmath.sym(0.5*(L+L.T))
  w_inc = nemlmath.skew(0.5*(L-L.T))
  
  s_n = np.zeros((6,))
  d_n = np.zeros((6,))
  w_n = np.zeros((3,))

  u_n = 0.0
  p_n = 0.0

  for i in range(steps):
    print(i)
    d_np1 = d_n + d_inc * dt
    w_np1 = w_n + w_inc * dt
    s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = pmodel.update_ld_inc(
        d_np1, d_n, w_np1, w_n, T, T, dt, 0,
        s_n, h_n, u_n, p_n)
    
    d_n = np.copy(d_np1)
    w_n = np.copy(w_np1)

    s_n = np.copy(s_np1)
    h_n = np.copy(h_np1)

    u_n = u_np1
    p_n = p_np1
  
  polefigures.pole_figure_discrete(orientations,[1,1,1],lattice)
  plt.title("Initial, <111>")
  plt.show()

  polefigures.pole_figure_discrete(pmodel.orientations(h_np1),[1,1,1],lattice)
  plt.title("Final, <111>")
  plt.show()

  polefigures.pole_figure_discrete(orientations,[1,0,0],lattice)
  plt.title("Initial, <100>")
  plt.show()

  polefigures.pole_figure_discrete(pmodel.orientations(h_np1),[1,0,0],lattice)
  plt.title("Final, <100>")
  plt.show()
