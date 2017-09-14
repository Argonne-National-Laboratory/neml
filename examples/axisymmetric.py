#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import axisym, neml, elasticity, surfaces

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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

  epoints = amodel.ls/2 + amodel.mesh[:-1]

  amodel.step(tup)

  stresses = np.array(amodel.stresses)
  plt.plot(epoints, np.mean(stresses[-1], axis=1)[:,0], 'k-')
  plt.plot(epoints, np.mean(stresses[-1], axis=1)[:,1], 'r-')
  plt.plot(epoints, np.mean(stresses[-1], axis=1)[:,2], 'b-')
  plt.show()

  strains = np.array(amodel.strains)
  plt.plot(epoints, np.mean(strains[-1], axis=1)[:,0], 'k-')
  plt.plot(epoints, np.mean(strains[-1], axis=1)[:,1], 'r-')
  plt.plot(epoints, np.mean(strains[-1], axis=1)[:,2], 'b-')
  plt.show()
