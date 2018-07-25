#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import axisym, neml, elasticity, surfaces, creep

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import scipy.interpolate as inter

def progress(f, l = 20):
  """
    Progress bar

    Parameters:
      f         completed fraction

    Optional:
      l         length
  """
  fill = int(round(f*l))
  pprint = "\rProgress: [%s] %6.2f%%" % ("#"*fill + "-"*(l-fill), f*100)
  sys.stdout.write(pprint)
  sys.stdout.flush()

def gen_material(E, nu, alpha):
  youngs = elasticity.YoungsModulus(E)
  poissons = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(youngs, poissons)

  return neml.SmallStrainElasticity(elastic, alpha = alpha)

if __name__ == "__main__":
  E = 100000.0
  nu = 0.3
  cte = 1.0e-3
  base = gen_material(E, nu, cte)

  Ro = 1750.0
  t = 50.0

  T1 = 0.0
  T2 = 100.0

  time = 100.0

  def temperature(r, t):
    return T1 + (T2-T1)/time * t

  def pressure(t):
    return 0.0
  
  times = np.linspace(0, time, 51)

  amodel = axisym.AxisymmetricProblem([Ro-t,Ro], [base], 
      [25], temperature, pressure, constrained = True)

  for t in times[1:]:
    amodel.step(t)
  
  stresses = np.array(amodel.stresses)
  print(np.mean(stresses[-1,:,0,2]))
  print(-(T2-T1) * E * cte) 
