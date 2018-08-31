#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, neml, elasticity, surfaces, hardening, ri_flow, uniaxial, axisym

import matplotlib.pyplot as plt
import numpy as np

def single_material_generator(E, nu, sY, alpha):
  youngs = elasticity.YoungsModulus(E)
  poisson = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(youngs, poisson)
  surface = surfaces.IsoJ2()
  model = neml.SmallStrainPerfectPlasticity(
    elastic, surface, sY, alpha = alpha)

  return model

if __name__ == "__main__":
  E = 100000.0
  nu = 0.25
  sY = 10000.0
  alpha = 1.0e-5
  n = 100

  mat = single_material_generator(E, nu, sY, alpha)

  pi = 10.0
  po = 5.0
  
  Ri = 100.0
  Ro = 120.0

  p_int = lambda t: t * pi
  p_ext = lambda t: t * po

  bree = axisym.BreeProblem([Ri, Ro], [mat], [n], lambda r, t: 0.0, p_int, p_ext = p_ext)

  bree.step(1.0)

  print("Bree problem:")
  print("'Hoop' stress")
  print("\tModel\t\t%f" % bree.stresses[-1][0])
  print("\tAnalytic\t%f" % ((pi * Ri - po * Ro) / (Ro - Ri)))
  print("\n")

  axi = axisym.AxisymmetricProblem([Ri,Ro], [mat], [n], lambda r, t: 0.0, p_int, 
      p_ext = p_ext)

  axi.step(1.0)
  
  s_inner = axi.stresses[-1][0][0]
  s_outer = axi.stresses[-1][-1][0]
  
  s_theta = lambda r: ((pi * Ri**2.0) / (Ro**2.0 - Ri**2.0) * (1.0 + Ro**2.0 / r**2.0) - 
      (po*Ro**2.0)/(Ro**2.0 - Ri**2.0) * (1.0 + Ri**2.0 / r**2.0))
  s_radial = lambda r: (pi * Ri**2.0 / (Ro**2.0 - Ri**2.0) * (1.0 - Ro**2.0 / r**2.0) -
      po * Ro**2.0 / (Ro**2.0 - Ri**2.0) * (1.0 - Ri**2.0 / r**2.0))
  s_axial = lambda r: pi * Ri**2.0 / (Ro**2.0 - Ri**2.0) - po * Ro**2.0 / (Ro**2.0 - Ri**2.0)
  
  print("Axisymmetric problem:")
  print("Inner")
  print("\tRadial stress")
  print("\t\tModel\t\t%f" % s_inner[0])
  print("\t\tAnalytic\t%f" % s_radial(Ri))

  print("\tHoop stress")
  print("\t\tModel\t\t%f" % s_inner[1])
  print("\t\tAnalytic\t%f" % s_theta(Ri))

  print("\tAxial stress")
  print("\t\tModel\t\t%f" % s_inner[2])
  print("\t\tAnalytic\t%f" % s_axial(Ri))

  print("Outer")
  print("\tRadial stress")
  print("\t\tModel\t\t%f" % s_outer[0])
  print("\t\tAnalytic\t%f" % -po)

  print("\tHoop stress")
  print("\t\tModel\t\t%f" % s_outer[1])
  print("\t\tAnalytic\t%f" % s_theta(Ro))

  print("\tAxial stress")
  print("\t\tModel\t\t%f" % s_outer[2])
  print("\t\tAnalytic\t%f" % s_axial(Ro))
