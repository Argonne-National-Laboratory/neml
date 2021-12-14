#!/usr/bin/env python3

import sys
sys.path.append('../..')

import os, glob
import scipy.interpolate as inter

import numpy as np
from neml import models, interpolate, elasticity, history
from neml.cp import hucocks, crystallography, sliprules, slipharden, inelasticity, kinematics, singlecrystal, polycrystal, polefigures, postprocessors
from neml.math import rotations, tensors, nemlmath, matrix
from neml import drivers


import matplotlib.pyplot as plt
import numpy.linalg as la
import numpy.random as ra
import pandas as pd
import xarray as xr
import tqdm
import warnings
warnings.filterwarnings("ignore")


def simplify_model(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s,
            k1_1, k1_2, k1_3, X, 
            k2_1, k2_2, k2_3,
            T = 298.0, emax = 0.05, N = 1, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True, PTR = True):  
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
  tau0 = np.array([taus_1]*3+[taus_2]*3+[taus_3]*6+[taut_1]*6+[taut_2]*6)
  
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

  mu = np.ones((24,))
  mu_slip = 30000.0
  mu_twin = 25000.0
  mu[:12] = mu_slip
  mu[12:] = mu_twin


  k1 = np.array([k1_1]*3+[k1_2]*3+[k1_3]*6)

  k2 = np.array([k2_1]*3+[k2_2]*3+[k2_3]*6)
  
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
  if PTR:
    single_model = singlecrystal.SingleCrystalModel(kmodel, lattice, 
        postprocessors = [], verbose = False, linesearch = True,
        initial_rotation = rotations.Orientation(0,0,0,angle_type="degrees"),
        miter = 100, max_divide = 10)
  else:
    single_model = singlecrystal.SingleCrystalModel(kmodel, lattice, 
        verbose = False, linesearch = True,
        initial_rotation = rotations.Orientation(0,0,0,angle_type="degrees"),
        miter = 100, max_divide = 10)
        
  if Taylor:
    orientations = rotations.random_orientations(N)
    model = polycrystal.TaylorModel(single_model, orientations, nthreads = nthreads)
    return drivers.uniaxial_test(model, strain_rate, T = T, 
            emax = emax, verbose = verbose)
  else:
    return drivers.uniaxial_test(single_model, strain_rate, 
            T = T, emax = emax, verbose = verbose)





def make_model(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s, 
            k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = 0.05, N = 1, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True,
            PTR = True):

  # unit transformer
  ut = 1.0e9
  
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
  # tau0 = np.array([170.0]*3+[90.5]*3+[210]*6+[180.0]*6+[250.0]*6)
  tau0 = np.array([taus_1]*3+[taus_2]*3+[taus_3]*6+[taut_1]*6+[taut_2]*6)
  
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

  mu = np.ones((24,))
  mu_slip = 30000.0
  mu_twin = 25000.0
  mu[:12] = mu_slip
  mu[12:] = mu_twin


  k1 = np.array([k1_1]*3+[k1_2]*3+[k1_3]*6)
  b = np.array([2.9511e-10]*3+[2.9511e-10]*3+[5.5364e-10]*6)*ut 
  gamma_dot = 1.0e7
  g = np.array([g_1]*3+[g_2]*3+[g_3]*6)
  tau_D = np.array([tau_D1]*3+[tau_D2]*3+[tau_D3]*6)
  eps_dot = strain_rate
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
  if PTR:
    single_model = singlecrystal.SingleCrystalModel(kmodel, lattice, 
        postprocessors = [], verbose = False, linesearch = True,
        initial_rotation = rotations.Orientation(0,0,0,angle_type="degrees"),
        miter = 100, max_divide = 10)
  else:
    single_model = singlecrystal.SingleCrystalModel(kmodel, lattice, 
        verbose = False, linesearch = True,
        initial_rotation = rotations.Orientation(0,0,0,angle_type="degrees"),
        miter = 100, max_divide = 10)
      
  if Taylor:
    orientations = rotations.random_orientations(N)
    model = polycrystal.TaylorModel(single_model, orientations, nthreads = nthreads)
    return drivers.uniaxial_test(model, strain_rate, T = T, 
            emax = emax, verbose = verbose)
  else:
    return drivers.uniaxial_test(single_model, strain_rate, 
            T = T, emax = emax, verbose = verbose)


def interpolate(strain, stress, targets):
  """
    This is just to make sure all our values line up, in case the model
    adaptively integrated or something
  """
  return inter.interp1d(strain, stress)(targets) 


def load_file(path):
  for file in glob.glob(path + "CG_Ti.csv"):
    df = pd.read_csv(file, usecols=[0,1], names=['Nominal_strain', 'True_stress'], header=None)
  return df



if __name__ == "__main__":
  
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  
  N = 1
  nthreads = 1
  strain_rate = 1.0e-4
  T = 298.0
  X_s = 0.9
  k1_1, k1_2, k1_3 = 1.50, 0.25, 5.0
  X = 0.3
  g_1, g_2, g_3 = 0.002, 0.002, 0.0055
  tau_D1, tau_D2, tau_D3 = 100.0, 100.0, 100.0
  
  res = make_model(X_s, k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, N = 1, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True)
            
            
  df = load_file(path_1)
  
  stress = interpolate(df['Nominal_strain'], df['True_stress'], res['strain'])
  
  plt.plot(res['strain'], res['stress'], label = "Sim - %3.0f C" % (T-273.15))
  plt.plot(res['strain'], stress, label = "Exp - %3.0f C" % (T-273.15))

  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  # plt.savefig("tension-Ti.png")
  plt.show()
  plt.close()




