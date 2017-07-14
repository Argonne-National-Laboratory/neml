#!/usr/bin/env python

import sys
sys.path.append('..')

import numpy as np

import matplotlib.pyplot as plt

from neml import solvers, neml, elasticity, drivers, surfaces, hardening, ri_flow, visco_flow, general_flow

def cycle_generator(de, edot):
  period = np.abs(2.0 * de/ edot)

  edot_sign = np.sign(de) * edot
  
  if de == 0.0:
    ediff = lambda t: t * 0.0 + T0
  else:
    ediff = lambda t: np.piecewise(t, [t < period / 2, t >= period / 2],
        [lambda tt: edot_sign * tt, lambda tt: de - edot_sign * (tt - period / 2)]) 

  return ediff, period

def grid_classical(max_x = 1.1, max_y = 2.5, A = 1.0, sigY = 100.0, 
    E = 10000.0, nu = 0.3, nx = 11, ny = 25, edot = 1.0e-3, 
    ncycles = 10, load_time = 1.0):
  """
    Simulate a classical two bar problem.
  """
  osx = max_x / (nx)
  osy = max_y / (ny)

  E_model = elasticity.YoungsModulus(E)
  nu_model = elasticity.PoissonsRatio(nu)
  emodel = elasticity.IsotropicLinearElasticModel(E_model, nu_model)

  surface = surfaces.IsoJ2()

  model = neml.SmallStrainPerfectPlasticity(emodel, surface, sigY)

  points = []
  conditions = []
  for x in np.linspace(osx/2, max_x-osx/2, nx):
    for y in np.linspace(osy/2, max_y-osy/2, ny):
      P = x * A * sigY * 2.0
      de = y * sigY / E * 2.0

      dstrain, period = cycle_generator(de, edot)
      res = drivers.twobar_test(model, A, A, lambda t: 0.0, lambda t: 0.0, 
          period, P, load_time, ncycles, dstrain = dstrain)

      points.append([x,y])
      conditions.append(res['classification'])

  return points, conditions, osx, osy, nx, ny

def plot_diagram(dx, dy, nx, ny, conditions):
  """
    Makes a pretty plot of the results
  """
  cmap = {
      "collapse": (1.0,0.0,0.0),
      "elastic" : (0.0,0.0,1.0),
      "elastic shakedown": (0.0,1.0,0.0),
      "plastic shakedown": (1.0,1.0,0.0),
      "ratcheting": (0.0,1.0,1.0)
      }
  img = np.zeros((ny, nx, 3))
  k = 0
  for i in range(nx):
    for j in range(ny):
      img[j,i] = cmap[conditions[k]]
      k += 1

  img = img[::-1,:,:]

  plt.imshow(img, interpolation = 'none')
  
  xmap = lambda x: x / dx - 0.5
  xpoints = [0.0, 0.5, 1]
  xpixels = [xmap(x) for x in xpoints]
  
  plt.xticks(xpixels, xpoints)

  ymap = lambda y: ny - 1 - (y / dy - 0.5)
  ypoints = [0.0, 1.0, 2.0]
  ypixels = [ymap(y) for y in ypoints]

  plt.yticks(ypixels, ypoints)

  plt.show()

if __name__ == "__main__":
  points, conditions, dx, dy, nx, ny = grid_classical()

  plot_diagram(dx, dy, nx, ny, conditions)

