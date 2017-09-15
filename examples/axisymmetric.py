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

  amodel = axisym.AxisymmetricProblem([240.0,242.0,250.0], 
      [gen_material(100000, 0.3, 1000.0, 20.0e-6),
        gen_material(100000, 0.3, 1000.0, 20.0e-6)], [25,25], 
      gradient, lambda t: p/tup * t, bias = True)

  amodel.step(tup)
  stresses = np.array(amodel.stresses)
  strains = np.array(amodel.strains)

  stress_vec = [stresses[-1,:,:,i].flatten() for i in range(6)]
  strain_vec = [strains[-1,:,:,i].flatten() for i in range(6)]

  epoints = np.array(amodel.ri).flatten()

  plt.plot(epoints, stress_vec[0], 'k-')
  plt.plot(epoints, stress_vec[1], 'r-')
  plt.plot(epoints, stress_vec[2], 'b-')
  plt.show()

  plt.plot(epoints, strain_vec[0], 'k-')
  plt.plot(epoints, strain_vec[1], 'r-')
  plt.plot(epoints, strain_vec[2], 'b-')
  plt.show()
