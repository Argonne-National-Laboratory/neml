#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures, postprocessors
from neml.math import rotations, tensors, nemlmath, matrix
from neml import elasticity, drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
  # Number of crystals and number of threads
  N = 1000
  nthreads = 3
  
  # Strain direction, rate, and number of steps
  L = np.array([[0.0,0,0],[0,1.0,0],[0,0,-1.0]])
  erate = 1.0e-4
  steps = 100
  emax = 2.0
  
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
  
  print(M)
  
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
  
  # Sets up the linear elastic tensor
  E = 100000.0
  nu = 0.3
  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  
  # emodel = elasticity.TransverseIsotropicLinearElasticModel(
      # C11,C33,C12,C13,C44,"components")
  
  # Sets up the slip system strength model (this is what you'll change)
  strength = slipharden.SimpleLinearHardening(M, tau0)
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
  

  res = drivers.uniaxial_test(tmodel, erate = erate, emax = emax, T = 298.0, verbose = True)
  true_strain = np.log(1+np.abs(res['strain']))
  true_stress = np.array(res['stress']) * (1 + np.abs(res['strain']))     
  plt.plot(true_strain, true_stress, 'g--', label = 'stress-strain')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.legend()
  plt.grid(True)
  plt.savefig("Dislocation-history-{}.png".format(int(T)))
  plt.show()
  plt.close()    