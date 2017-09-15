#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import axisym, neml, elasticity, surfaces

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import scipy.interpolate as inter

def gen_material(E, nu, sY, alpha):
  """
    Generate a perfectly plastic material model.
  """
  youngs = elasticity.YoungsModulus(E)
  poissons = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(youngs, poissons)
  surface = surfaces.IsoJ2()

  return neml.SmallStrainPerfectPlasticity(elastic, surface, sY, alpha = alpha)

if __name__ == "__main__":
  tup = 10.0
  gradient = axisym.generate_thickness_gradient(240.0, 250.0, 400, 600,
      10.0, hold = 100.0, Tdot_cold = 20.0, hold_together = 20.0,
      delay = tup)
  p = 1.1
  period = 150

  def pressure(t):
    if t < tup:
      return p/tup * t
    else:
      return p

  amodel = axisym.AxisymmetricProblem([240.0,242.0,250.0], 
      [gen_material(100000, 0.3, 150.0, 20.0e-6),
        gen_material(100000, 0.3, 150.0, 20.0e-6)], [50,50], 
      gradient, pressure, bias = False)
  
  ncycles = 5
  total = tup + period * ncycles

  for t in np.linspace(0,total, 500):
    amodel.step(t, verbose = False)

  """
  stresses = np.array(amodel.stresses)
  strains = np.array(amodel.strains)
  tstrains = np.array(amodel.tstrains)

  stress_vec = [stresses[-1,:,:,i].flatten() for i in range(6)]
  strain_vec = [strains[-1,:,:,i].flatten() for i in range(6)]
  tstrain_vec = [tstrains[-1,:,:,i].flatten() for i in range(6)]

  epoints = np.array(amodel.ri).flatten()

  plt.plot(epoints, stress_vec[0], 'k-')
  plt.plot(epoints, stress_vec[1], 'r-')
  plt.plot(epoints, stress_vec[2], 'b-')
  plt.show()

  plt.plot(epoints, strain_vec[0], 'k-')
  plt.plot(epoints, strain_vec[1], 'r-')
  plt.plot(epoints, strain_vec[2], 'b-')
  plt.show()

  plt.plot(epoints, tstrain_vec[0], 'k-')
  plt.plot(epoints, tstrain_vec[1], 'r-')
  plt.plot(epoints, tstrain_vec[2], 'b-')
  plt.show()
  """
