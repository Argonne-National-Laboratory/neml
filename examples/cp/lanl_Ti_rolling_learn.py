#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures, postprocessors
from neml.math import rotations, tensors, nemlmath, matrix
from neml import elasticity
from neml import drivers

import numpy.linalg as la
import numpy.random as ra

import tqdm
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
  
  # unit transformer
  ut = 1.0e9

  # Temperature in K
  T = 298.0
  
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
  
  # Reference slip rate and rate sensitivity exponent
  g0 = 1.0
  n = 12.0
  
  # Twin threshold 
  twin_threshold = 0.75
  
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
  
  # sets up input parameters
  
  # Sets up the interaction matrix
  num_basal, num_prism, num_pyram = 3, 3, 6
  num_ttwin, num_ctwin = 6, 6
  ## basal plane
  C1_single = np.array([50.0]*num_ttwin+[50.0]*num_ttwin)
  C1 = np.stack([C1_single for _ in range(num_basal)])
  ## prismatic plane
  C2_single = np.array([100.0]*num_ttwin+[1000.0]*num_ttwin)
  C2 = np.stack([C2_single for _ in range(num_prism)])
  ## pyramidal plane
  C3_single = np.array([250.0]*num_ttwin+[170.0]*num_ttwin)
  C3 = np.stack([C3_single for _ in range(num_pyram)])
  ## stack up each slip system
  C_np = np.vstack((C1, C2, C3)).T

  C_st = matrix.SquareMatrix(12, type = "dense", data = C_np.flatten())

  M = matrix.SquareMatrix(24, type = "diagonal_blocks",
      data = [10.0, 10.0], blocks = [12,12])

  mu = np.ones((24,))
  mu_slip = 43000.0 #30000.0
  mu_twin = 43000.0 #25000.0
  mu[:12] = mu_slip
  mu[12:] = mu_twin
  X_s = 0.9
  k1 = np.array([1.50e9]*3+[2.5e8]*3+[5.00e9]*6)/ut #+[1.50]*3+[0.25]*3+[5.00]*6)
  X = 0.3
  b = np.array([2.9511e-10]*3+[2.9511e-10]*3+[5.5364e-10]*6)*ut  #+[0.29511]*3+[0.29511]*3+[0.55364]*6)
  # b = np.array([0.29511]*3+[0.29511]*3+[0.55364]*6+[0.060040]*6+[0.082015]*6)
  gamma_dot = 1.0e7
  g = np.array([0.002]*3+[0.002]*3+[0.0055]*6) #+[0.002]*3+[0.002]*3+[0.0055]*6)
  tau_D = np.array([100.0]*3+[100.0]*3+[100.0]*6) #+[100.0]*3+[100.0]*3+[100.0]*6)
  eps_dot = erate
  k = 1.38064852e-23*ut**2.0

  k2 = k1*X*b*(1-k*T/(tau_D*b**3)*np.log(eps_dot/gamma_dot))/g

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
      postprocessors = [], verbose = True, linesearch = True,
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

  for i in tqdm.tqdm(range(steps)):
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
  plt.close()


  res = drivers.uniaxial_test(pmodel, erate, T = T, verbose = True) 
  plt.plot(res['strain'], res['stress'], label = "%3.0f C" % (T-273.15))

  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  # plt.savefig("tension-Ti.png")
  plt.show()
  plt.close()