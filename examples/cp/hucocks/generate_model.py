#!/usr/bin/env python3

import numpy as np

import sys
sys.path.append('../../..')

from neml import models, interpolate, elasticity, history
from neml.cp import hucocks, crystallography, sliprules, slipharden, inelasticity, kinematics, singlecrystal, polycrystal
from neml.math import rotations

def make_model(N, nthreads = 1):
  Ts = np.array([500.0,550.0,600.0,650.0]) + 273.15

  L = crystallography.CubicLattice(1.0)
  L.add_slip_system([1,1,0],[1,1,1])

  uf = 1.0e-9

  J1 = 2e14 * uf**2.0
  J2 = 3.3e14 * uf**2.0
  K = 2.56e-30 / uf**4.0
  L0 = 3.1623e-7 / uf
  b = 2.5e-10
  b_d = b / uf
  ad = 0.35
  G = interpolate.PiecewiseLinearInterpolate(
      list(Ts),
      [61068, 59541.0, 57633.6, 55725.2])
  
  dmodel = hucocks.DislocationSpacingHardening(J1, J2, K, L0, ad, 
      b_d, G, L)

  # Setup for [Cr,C] <-> M23C6
  am_car = 3.6e-10
  N0_car = 1.0e13
  Vm_car = 6e-6
  chi_car = 0.3
  D0_car = 1.5e-4
  Q0_car = 240e3
  c0_car = [interpolate.ConstantInterpolate(16.25/100.0),
      interpolate.ConstantInterpolate(0.0375/100.0)]
  cp_car = [interpolate.PiecewiseLinearInterpolate(list(Ts),
    [69.85/100, 69.05/100, 68.32/100, 67.52/100]),
    interpolate.PiecewiseLinearInterpolate(list(Ts),
      [5.13/100, 5.13/100, 5.13/100, 5.13/100])]
  ceq_car = [interpolate.PiecewiseLinearInterpolate(list(Ts),
    [15.64/100,15.69/100,15.75/100,15.83/100]),
    interpolate.PiecewiseLinearInterpolate(list(Ts),
      [7.25e-6/100, 2.92e-5/100, 9.48e-5/100, 2.97e-4/100])]
  Cf_car = interpolate.PiecewiseLinearInterpolate(list(Ts), 
      [1.0, 1.0, 0.3, 0.03])

  carbide = hucocks.HuCocksPrecipitationModel(c0_car, cp_car, ceq_car, 
      am_car, N0_car, Vm_car, chi_car, D0_car,
      Q0_car, Cf_car) 

  am_laves = 3.6e-10
  N0_laves = 5e14
  Vm_laves = 2e-6
  chi_laves = 0.25
  D0_laves = 7.4e-4
  Q0_laves = 283e3
  c0_laves = [2.33/100.0]
  cp_laves = [50.0/100.0]
  ceq_laves = [interpolate.PiecewiseLinearInterpolate(list(Ts),
    [0.25/100,0.46/100.0,0.76/100.0,1.16/100.0])]
  Cf_laves = 1.0

  laves = hucocks.HuCocksPrecipitationModel(c0_laves, cp_laves, ceq_laves, 
      am_laves, N0_laves, Vm_laves, chi_laves, D0_laves,
      Q0_laves, Cf_laves) 

  ap = 0.84
  ac = 0.000457

  tau_model = hucocks.HuCocksHardening(dmodel, [carbide, laves], ap, ac, b, G)

  g0 = 1.0
  a0 = 0.5
  G0 = 77000.0e6
  A = 3.0/4.0
  B = 4.0/3.0

  slip_model = hucocks.ArrheniusSlipRule(tau_model, g0, A, B, b,
      a0, G0)

  imodel = inelasticity.AsaroInelasticity(slip_model)
  
  youngs = interpolate.PiecewiseLinearInterpolate(list(Ts),
      [160000.0, 156000.0, 151000.0, 140000.0])
  emodel = elasticity.IsotropicLinearElasticModel(youngs, "youngs", 
      0.31, "poissons")
  kmodel = kinematics.StandardKinematicModel(emodel, imodel)
  smodel = singlecrystal.SingleCrystalModel(kmodel, L, verbose = True)

  orientations = rotations.random_orientations(N)

  model = polycrystal.TaylorModel(smodel, orientations, nthreads = nthreads)

  return smodel

if __name__ == "__main__":
  model = make_model(100)

  
