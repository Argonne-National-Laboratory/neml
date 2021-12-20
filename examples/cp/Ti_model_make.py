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


def Ti_maker_sim(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s, 
            k1_1, k1_2, k1_3,
            k2_1, k2_2, k2_3,
            T = 296.0, emax = 0.05, N = 1, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True,
            PTR = True):

  # temperature levels
  Ts = np.array([296.0, 373.0, 473.0, 573.0, 673.0, 773.0])
  # unit transformer
  ut = 1.0e9
  
  # Model
  a = 2.9511*0.1 # nm
  c = 4.68433*0.1 # nm
  
  # Elastic constants in MPa
  C11 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [160000.0, 157900.0, 152200.0, 146800.0, 141600.0, 136800.0])
  C33 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [181000.0, 177400.0, 173400.0, 169600.0, 166100.0, 162700.0])
  C44 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [46500.0, 45300.0, 43400.0, 41400.0, 39200.0, 37000.0])      
  C12 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [90000.0, 93400.0, 95200.0, 96700.0, 97800.0, 98500.0])
  C13 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [66000.0, 69400.0, 69500.0, 69200.0, 69000.0, 68800.0])
  
  
  # Constant part of the strength for slip and twin
  taus_1 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_1, taus_1, taus_1, taus_1, taus_1, taus_1]) 
  taus_2 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_2, taus_2, taus_2, taus_2, taus_2, taus_2]) 
  taus_3 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_3, taus_3, taus_3, taus_3, taus_3, taus_3]) 
  taut_1 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taut_1, taut_1, taut_1, taut_1, taut_1, taut_1]) 
  taut_2 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taut_2, taut_2, taut_2, taut_2, taut_2, taut_2]) 
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


  # source from Fracture of Titanium Alloys at High Strain Rates and under Stress Triaxiality
  # calculate temperature depdendent shear modulus of slip systems  u = 39.61-0.03223*T
  mu_slip = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [30005.46, 27588.21, 24365.21, 21142.21, 17919.21, 14696.21])
  # calculate temperature depdendent shear modulus of twin systems  u = 34.605-0.03223*T
  mu_twin = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [25000.46, 22583.21, 19360.21, 16137.21, 12914.21, 9691.21])
  
  mu = np.array([mu_slip]*12+[mu_twin]*12)
  
  k1 = np.array([k1_1]*3+[k1_2]*3+[k1_3]*6)  # k1_1, k1_2, k1_3 = 1.35, 0.25, 5.00

  k2_1 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [k2_1, k2_1, k2_1, k2_1, k2_1, k2_1]) 
  k2_2 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [k2_2, k2_2, k2_2, k2_2, k2_2, k2_2]) 
  k2_3 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [k2_3, k2_3, k2_3, k2_3, k2_3, k2_3]) 

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
            emax = emax, sdir = np.array([-1,0,0,0,0,0]), verbose = verbose)
  else:
    return drivers.uniaxial_test(single_model, strain_rate, 
            T = T, emax = emax, sdir = np.array([-1,0,0,0,0,0]), verbose = verbose)

def Ti_maker(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s, 
            k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = 0.05, N = 1, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True,
            PTR = True):

  # temperature levels
  Ts = np.array([298.0, 373.0, 473.0, 573.0, 673.0, 773.0])
  # unit transformer
  ut = 1.0e9
  
  # Model
  a = 2.9511*0.1 # nm
  c = 4.68433*0.1 # nm
  
  # Elastic constants in MPa
  C11 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [160000.0, 157900.0, 152200.0, 146800.0, 141600.0, 136800.0])
  C33 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [181000.0, 177400.0, 173400.0, 169600.0, 166100.0, 162700.0])
  C44 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [46500.0, 45300.0, 43400.0, 41400.0, 39200.0, 37000.0])      
  C12 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [90000.0, 93400.0, 95200.0, 96700.0, 97800.0, 98500.0])
  C13 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [66000.0, 69400.0, 69500.0, 69200.0, 69000.0, 68800.0])
  
  
  # Constant part of the strength for slip and twin
  taus_1 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_1, taus_1, taus_1, taus_1, taus_1, taus_1]) 
  taus_2 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_2, taus_2, taus_2, taus_2, taus_2, taus_2]) 
  taus_3 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_3, taus_3, taus_3, taus_3, taus_3, taus_3]) 
  taut_1 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taut_1, taut_1, taut_1, taut_1, taut_1, taut_1]) 
  taut_2 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taut_2, taut_2, taut_2, taut_2, taut_2, taut_2]) 
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


  # source from Fracture of Titanium Alloys at High Strain Rates and under Stress Triaxiality
  # calculate temperature depdendent shear modulus of slip systems  u = 39.61-0.03223*T
  mu_slip = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [30005.46, 27588.21, 24365.21, 21142.21, 17919.21, 14696.21])
  # calculate temperature depdendent shear modulus of twin systems  u = 34.605-0.03223*T
  mu_twin = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [25000.46, 22583.21, 19360.21, 16137.21, 12914.21, 9691.21])
  
  mu = np.array([mu_slip]*12+[mu_twin]*12)
  
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
        # initial_rotation = rotations.Orientation(0,0,0,angle_type="degrees"),
        miter = 100, max_divide = 10)
  else:
    single_model = singlecrystal.SingleCrystalModel(kmodel, lattice, 
        verbose = False, linesearch = True,
        # initial_rotation = rotations.Orientation(0,0,0,angle_type="degrees"),
        miter = 100, max_divide = 10)
      
  if Taylor:
    orientations = rotations.random_orientations(N)
    model = polycrystal.TaylorModel(single_model, orientations, nthreads = nthreads)
    return drivers.uniaxial_test(model, strain_rate, T = T, 
            emax = emax, verbose = verbose)
  else:
    return drivers.uniaxial_test(single_model, strain_rate, 
            T = T, emax = emax, verbose = verbose)


def Ti_orientation(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s, 
            k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = 0.05, N = 1, 
            strain_rate = 1.0e-4, nthreads = 1, 
            steps = 100, verbose = True, Taylor = True,
            PTR = True):

  # temperature levels
  Ts = np.array([298.0, 373.0, 473.0, 573.0, 673.0, 773.0])
  # unit transformer
  ut = 1.0e9
  
  # Model
  a = 2.9511*0.1 # nm
  c = 4.68433*0.1 # nm
  
  # Elastic constants in MPa
  C11 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [160000.0, 157900.0, 152200.0, 146800.0, 141600.0, 136800.0])
  C33 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [181000.0, 177400.0, 173400.0, 169600.0, 166100.0, 162700.0])
  C44 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [46500.0, 45300.0, 43400.0, 41400.0, 39200.0, 37000.0])      
  C12 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [90000.0, 93400.0, 95200.0, 96700.0, 97800.0, 98500.0])
  C13 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [66000.0, 69400.0, 69500.0, 69200.0, 69000.0, 68800.0])
  
  
  # Constant part of the strength for slip and twin
  taus_1 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_1, taus_1, taus_1, taus_1, taus_1, taus_1]) 
  taus_2 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_2, taus_2, taus_2, taus_2, taus_2, taus_2]) 
  taus_3 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taus_3, taus_3, taus_3, taus_3, taus_3, taus_3]) 
  taut_1 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taut_1, taut_1, taut_1, taut_1, taut_1, taut_1]) 
  taut_2 = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [taut_2, taut_2, taut_2, taut_2, taut_2, taut_2]) 
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
  
  # Sets up the actual deformation tensor
  L = np.array([[1.0,0,0],[0,0.0,0],[0,0,0.0]])
  L *= erate
  dt = emax / steps / erate
  # Randomly selects initial orientations
  orientations = rotations.random_orientations(N)
  
  # Plots an initial basal pole figure
  polefigures.pole_figure_discrete(orientations,[0,0,0,1],lattice)
  plt.title("Initial, <0001>")
  plt.savefig("Ti-initial-orientation.png")
  plt.show()
  plt.close()
  
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


  # source from Fracture of Titanium Alloys at High Strain Rates and under Stress Triaxiality
  # calculate temperature depdendent shear modulus of slip systems  u = 39.61-0.03223*T
  mu_slip = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [30005.46, 27588.21, 24365.21, 21142.21, 17919.21, 14696.21])
  # calculate temperature depdendent shear modulus of twin systems  u = 34.605-0.03223*T
  mu_twin = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [25000.46, 22583.21, 19360.21, 16137.21, 12914.21, 9691.21])
  
  mu = np.array([mu_slip]*12+[mu_twin]*12)
  
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
    pmodel = polycrystal.TaylorModel(single_model, orientations, nthreads = nthreads)
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
    plt.savefig("Ti-deformed-orientation.png")
    return plt.show(), plt.close() 
  else:
    return None
  


def interp(strain, stress, targets):
  """
    This is just to make sure all our values line up, in case the model
    adaptively integrated or something
  """
  return inter.interp1d(strain, stress)(targets) 


def load_file(path, temper):
  for file in glob.glob(path + temper + "k.csv"):
    df = pd.read_csv(file, usecols=[0,1], names=['Nominal_strain', 'True_stress'], header=None)
  return df

# def load_file(path, temper):
  # for file in glob.glob(path + "CG_Ti.csv"):
    # df = pd.read_csv(file, usecols=[0,1], names=['Nominal_strain', 'True_stress'], header=None)
  # return df


if __name__ == "__main__":
  
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  Ts = [298.0]
  
  for T in Ts:
    N = 10
    nthreads = 1
    erate = 1.0e-4
    emax = 0.15

    taus_1, taus_2, taus_3 = 270.0, 190.0, 310.0
    taut_1, taut_2 = 280.0, 350.0
    X_s = 0.9
    k1_1, k1_2, k1_3 = 2.00, 2.00, 5.00
    k2_1, k2_2, k2_3 = 70.0, 70.0, 70.0

    res = Ti_maker_sim(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s, 
            k1_1, k1_2, k1_3,
            k2_1, k2_2, k2_3,
            T = T, emax = emax, N = N, 
            strain_rate = erate, nthreads = nthreads, 
            verbose = True, Taylor = True,
            PTR = True)
            
            
    df = load_file(path_1, str(int(T)))

    plt.plot(res['strain'], res['stress'], 'g--', label = "Sim - %3.0f C" % (T-273.0))
    plt.plot(df['Nominal_strain'], df['True_stress'], 'b-', label = "Exp - %3.0f C" % (T-273.0))
    plt.grid(True)
    plt.legend(loc='best')
    plt.xlabel("Strain (mm/mm)")
    plt.ylabel("Stress (MPa)")
    plt.savefig("tension-Ti-{}.png".format(T))
    # plt.show()
    plt.close()
  
  """
  for T in Ts:
    N = 1
    nthreads = 1
    erate = 1.0e-4
    emax = 0.15

    taus_1, taus_2, taus_3 = 270.0, 190.0, 310.0
    taut_1, taut_2 = 280.0, 350.0
    X_s = 0.9
    k1_1, k1_2, k1_3 = 2.00, 2.00, 5.00
    X = 0.35
    g_1, g_2, g_3 = 0.0045, 0.0045, 0.005
    tau_D1, tau_D2, tau_D3 = 100.0, 100.0, 100.0

    res = Ti_maker(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s, 
            k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = T, emax = emax, N = N, 
            strain_rate = erate, nthreads = nthreads, 
            verbose = True, Taylor = True,
            PTR = True)
            
            
    df = load_file(path_1, str(int(T)))

    plt.plot(res['strain'], res['stress'], 'g--', label = "Sim - %3.0f C" % (T-273.0))
    plt.plot(df['Nominal_strain'], df['True_stress'], 'b-', label = "Exp - %3.0f C" % (T-273.0))
    plt.grid(True)
    plt.legend(loc='best')
    plt.xlabel("Strain (mm/mm)")
    plt.ylabel("Stress (MPa)")
    plt.savefig("tension-Ti-{}.png".format(T))
    # plt.show()
    plt.close()
  """

  """
  # ) visualize the orientation evolution
  N = 10
  nthreads = 1
  erate = 1.0e-4
  emax = 2.0
  steps = 100
  T = 298.0

  taus_1, taus_2, taus_3 = 270.0, 190.0, 310.0
  taut_1, taut_2 = 280.0, 350.0
  X_s = 0.9
  k1_1, k1_2, k1_3 = 2.00, 2.00, 5.00
  X = 0.35
  g_1, g_2, g_3 = 0.0045, 0.0045, 0.005
  tau_D1, tau_D2, tau_D3 = 100.0, 100.0, 100.0

  Ti_orientation(taus_1, taus_2, taus_3,
        taut_1, taut_2, X_s, 
        k1_1, k1_2, k1_3, X, 
        g_1, g_2, g_3,
        tau_D1, tau_D2, tau_D3,
        T = T, emax = emax, N = N, 
        strain_rate = erate, nthreads = nthreads, 
        steps = steps, verbose = True, Taylor = True,
        PTR = True)
  """

