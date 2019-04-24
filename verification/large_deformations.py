#!/usr/bin/env python3

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

from neml import models, drivers, elasticity

def shear():
  """
    Recreate some of the results of Bazant 2013
  """
  F = lambda t: np.array([[1,t,0],[0,1.0,0],[0,0,1]])
  tmax = 10.0
  nsteps = 250
  
  G = 100.0
  K = 10.0

  elastic = elasticity.IsotropicLinearElasticModel(G, "shear", K, "bulk")
  model = models.SmallStrainElasticity(elastic)

  res = drivers.def_grad_driver(model, F, tmax, nsteps)
 
  theta = [Fi[0,1] / 2.0 for Fi in res['F']]
  nstress = [Si[5] / G / np.sqrt(2) for Si in res['stress']]
  plt.plot(theta, nstress, label = "Truesdell")

  model = models.SmallStrainElasticity(elastic, truesdell = False)
  
  res = drivers.def_grad_driver(model, F, tmax, nsteps)
  theta = [Fi[0,1] / 2.0 for Fi in res['F']]
  nstress = [Si[5] / G / np.sqrt(2) for Si in res['stress']]
  plt.plot(theta, nstress, label = "Jaumann")

  plt.xlabel(r"$\tan\left(\theta\right / 2)$")
  plt.ylabel(r"$\tau_{12}/G$")
  plt.legend(loc = 'best')
  plt.show()

def rotation():
  """
    Rotate a fixed stress tensor
  """
  def F(t):
    cut = 1.0
    bF = lambda tt: np.array([[1+tt,0,0],[0,1.0-tt/4.0,0],[0,0,1-tt/4.0]]) 
    R = lambda tt: np.array([
      [np.cos(tt*(np.pi/2)), -np.sin(tt*(np.pi/2)),0],
      [np.sin(tt*(np.pi/2)), np.cos(tt*(np.pi/2)),0],
      [0, 0, 1]])

    if t <= cut:
      F = bF(t)
    else:
      F = np.dot(R(t-cut), bF(cut))

    return F

  tmax = 2.0
  nsteps = 1000

  E = 100.0
  nu = 0.25

  elastic = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  model = models.SmallStrainElasticity(elastic)

  res = drivers.def_grad_driver(model, F, tmax, 2*nsteps)

  print(res['stress'][nsteps])
  print(res['stress'][2*nsteps])

if __name__ == "__main__":
  rotation()
  shear()
