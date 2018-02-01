#!/usr/bin/env python

import sys
sys.path.append('..')
from neml import neml, elasticity, arbbar

import numpy as np

def make_mat(E, a, nu = 0.3):
  youngs = elasticity.YoungsModulus(E)
  poissons = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(youngs, poissons)
  
  return neml.SmallStrainElasticity(elastic, alpha = a)

if __name__ == "__main__":
  T1 = 0
  T2 = 100.0
  E1 = 100000.0
  E2 = 150000.0
  E3 = 200000.0
  A1 = 1.5
  A2 = 1.25
  A3 = 1.0
  l1 = 11.0
  l2 = 5.0
  l3 = 8.0
  a1 = 1.0e-5
  a2 = 1.7e-5
  a3 = 1.9e-5

  expected = (T2-T1)*E1*E2*A1*A2*(a1*l1 - a2*l2 - a3*l3) / (
      A2*A3*l1*E2*E3 + A1*E1*(A2*l3*E2+A3*l2*E3))

  mat1 = make_mat(E1, a1)
  mat2 = make_mat(E2, a2)
  mat3 = make_mat(E3, a3)
  
  tload = 1.0
  T = lambda t: T1 + (T2-T1)/tload * t

  model = arbbar.BarModel()

  model.add_node(1)
  model.add_node(2)
  model.add_node(3)
  
  model.add_edge(1,3, object = arbbar.Bar(mat1, A1, l1, T = T))
  model.add_edge(1,2, object = arbbar.Bar(mat2, A2, l2, T = T))
  model.add_edge(2,3, object = arbbar.Bar(mat3, A3, l3, T = T))

  model.add_displacement_bc(3, lambda t: 0.0)

  nsteps = 1
  dt = tload / nsteps
  for i in range(nsteps):
    model.solve(dt)
  
  calced = model[2][3][0]['object'].mstrain[-1]

  print("Expected: %f" % expected)
  print("Calculated: %f" % calced)
