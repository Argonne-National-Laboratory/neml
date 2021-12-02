#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures, postprocessors
from neml.math import rotations, tensors, nemlmath, matrix
from neml import elasticity

import numpy.linalg as la
import numpy.random as ra

import matplotlib.pyplot as plt

if __name__ == "__main__":
  # Number of crystals and number of threads
  N = 100
  nthreads = 3
  
  # Strain direction, rate, and number of steps
  L = np.array([[0.0,0,0],[0,1.0,0],[0,0,-1.0]])
  erate = 1.0e-4
  steps = 100
  emax = 2.0
  
  # Temperature in K
  T = 0
  
  # Model
  a = 2.9511*0.1 # nm
  c = 4.68433*0.1 # nm
  
  # Elastic constants in MPa
  C11 = 160000.0
  C33 = 181000.0
  C44 = 46500.0
  C12 = 90000.0
  C13 = 66000.0
  
  # Constant part of the strength for slip and twin
  tau0 = np.array([170.0]*3+[90.5]*3+[210]*6+[180.0]*6+[250.0]*6)/10.0
  
  # Not realistic but I just want it to roll quickly
  # Hardening coefficients for slip (H1) and twinning (H2)
  H1 = 10.0
  H2 = 10.0
  
  # Reference slip rate and rate sensitivity exponent
  g0 = 1.0
  n = 12.0
  
  # Twin threshold 
  twin_threshold = 0.75
  
  # Sets up the interaction matrix
  M = matrix.SquareMatrix(24, type = "diagonal_blocks", 
      data = [H1,H2], blocks = [12,12])
  
  # Sets up the lattice crystallography
  lattice = crystallography.HCPLattice(a, c)
  # Basal <a>
  lattice.add_slip_system([1,1,-2,0],[0,0,0,1])
  # Prismatic <a>
  lattice.add_slip_system([1,1,-2,0],[1,0,-1,0])
  # Pyramidal <c+a>
  lattice.add_slip_system([1,1,-2,-3],[1,1,-2,2])
  # Tension twinning
  lattice.add_twin_system([-1,0,1,1],[1,0,-1,2],[1,0,-1,1],[1,0,-1,-2])
  # Compression twinning
  lattice.add_twin_system([1,1,-2,-3],[1,1,-2,2],[2,2,-4,3],[1,1,-2,-4])

  # Sets up the actual deformation tensor
  L *= erate
  dt = emax / steps / erate
  # Randomly selects initial orientations
  orientations = rotations.random_orientations(N)
  
  # Plots an initial basal pole figure
  polefigures.pole_figure_discrete(orientations,[0,0,0,1],lattice)
  plt.title("Initial, <0001>")
  plt.show()
  
  G_np = ra.random((12,12)) * 10
  C_st = matrix.SquareMatrix(12, type = "dense", data = G_np.flatten())
  print(C_st)
  mu = np.ones((24,))
  mu_slip = 30000.0
  mu_twin = 25000.0
  mu[:12] = mu_slip
  mu[12:] = mu_twin
  X_s = 0.9
  
  k1_v = 0.5
  k2_v = 0.75

  k1 = np.ones((12,)) * k1_v
  k2 = np.ones((12,)) * k2_v
  # Sets up the linear elastic tensor
  emodel = elasticity.TransverseIsotropicLinearElasticModel(
      C11,C33,C12,C13,C44,"components")
  
  # Sets up the slip system strength model (this is what you'll change)
  strength = slipharden.LANLTiModel(tau0, C_st, mu, k1, k2, X_s=X_s)
  # strength = slipharden.SimpleLinearHardening(M, tau0)
  # Sets up the slip rule
  slipmodel = sliprules.PowerLawSlipRule(strength, g0, n)
  # Sets up the model inelastic rate kinematics 
  imodel = inelasticity.AsaroInelasticity(slipmodel)
  # Sets up the overall model kinematics
  kmodel = kinematics.StandardKinematicModel(emodel, imodel)
  
  # This is the object that causes twins to recrystallize
  twinner = postprocessors.PTRTwinReorientation(twin_threshold)
  
  # Sets up the single crystal model
  model = singlecrystal.SingleCrystalModel(kmodel, lattice, 
      postprocessors = [twinner], verbose = False, linesearch = True,
      miter = 100, max_divide = 10)
  
  # Sets up the poly crystal model
  pmodel = polycrystal.TaylorModel(model, orientations, nthreads = nthreads)
  
  # Runs the exmaple
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
  
  # Plots a second, as-rolled basal pole figure
  polefigures.pole_figure_discrete(pmodel.orientations(h_np1),
      [0,0,0,1],lattice)
  plt.title("Final, <0001>")
  plt.show()
