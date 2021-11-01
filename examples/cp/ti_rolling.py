#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures, postprocessors
from neml.math import rotations, tensors, nemlmath, matrix
from neml import elasticity

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 200
  nthreads = 1

  L = np.array([[0,0,0],[0,1.0,0],[0,0,-1.0]])
  erate = 1.0e-4
  steps = 100
  emax = 1.0

  T = 0
  
  # Model
  a = 2.9511*0.1 # nm
  c = 4.68433*0.1 # nm

  C11 = 160000.0
  C33 = 181000.0
  C44 = 46500.0
  C12 = 90000.0
  C13 = 66000.0

  tau0 = np.array([170.0]*3+[90.5]*3+[210]*6+[180.0]*6+[250.0]*6)
  
  # Not realistic but I just want it to roll quickly
  H1 = 100.0 
  H2 = 100.0

  g0 = 1.0
  n = 12.0

  twin_threshold = 0.5

  M = matrix.SquareMatrix(24, type = "diagonal_blocks", 
      data = [H1,H2], blocks = [12,12])

  lattice = crystallography.HCPLattice(a, c)
  # Basal <a>
  lattice.add_slip_system([1,1,-2,0],[0,0,0,1])
  # Prismatic <a>
  lattice.add_slip_system([1,1,-2,0],[1,0,-1,0])
  # Pyramidal <c+a>
  lattice.add_slip_system([1,1,-2,3],[1,1,-2,2])
  # Tension twinning
  lattice.add_twin_system([-1,0,1,1],[1,0,-1,2],[1,0,-1,1],[1,0,-1,-2])
  # Compression twinning
  lattice.add_twin_system([1,1,-2,-3],[1,1,-2,2],[2,2,-4,3],[1,1,-2,-4])

  # Setup
  L *= erate
  dt = emax / steps / erate
  orientations = rotations.random_orientations(N)

  polefigures.pole_figure_discrete(orientations,[0,0,0,1],lattice)
  plt.title("Initial, <0001>")
  plt.show()

  emodel = elasticity.TransverseIsotropicLinearElasticModel(
      C11,C33,C12,C13,C44,"components")

  strength = slipharden.SimpleLinearHardening(M, tau0)
  slipmodel = sliprules.PowerLawSlipRule(strength, g0, n)
  imodel = inelasticity.AsaroInelasticity(slipmodel)
  kmodel = kinematics.StandardKinematicModel(emodel, imodel)

  twinner = postprocessors.PTRTwinReorientation(twin_threshold)

  model = singlecrystal.SingleCrystalModel(kmodel, lattice, 
      postprocessors = [twinner], verbose = True, linesearch = True)

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

  polefigures.pole_figure_discrete(pmodel.orientations(h_np1),
      [0,0,0,1],lattice)
  plt.title("Final, <0001>")
  plt.show()
