#!/usr/bin/env python

import sys
sys.path.append('..')

import numpy as np

import matplotlib.pyplot as plt

from neml import solvers, neml, elasticity, drivers, surfaces, hardening, ri_flow, visco_flow, general_flow

def cycle_generator(dT, T0, Tdot):
  period = np.abs(2.0 * dT/ Tdot)

  Tdot_sign = np.sign(dT) * Tdot
  
  T2 = lambda t: t * 0.0 + T0

  if dT == 0.0:
    T1 = lambda t: t * 0.0 + T0
  else:
    T1 = lambda t: np.piecewise(t, [t < period / 2, t >= period / 2],
        [lambda tt: T0 + Tdot_sign * tt, lambda tt: T0 + dT - Tdot_sign * (tt - period / 2)]) 

  return T1, T2, period

def grid_classical(g, max_x = 1.1, max_y = 2.5, A = 1.0, sigY = 100.0, 
    E = 10000.0, nu = 0.3, alpha = 1.0e-4, nx = 11, ny = 25, T0 = 0.0, Tdot = 1.0,
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

  model = neml.SmallStrainPerfectPlasticity(emodel, surface, sigY, alpha = alpha)

  points = []
  conditions = []
  for x in np.linspace(osx/2, max_x-osx/2, nx):
    for y in np.linspace(osy/2, max_y-osy/2, ny):
      P = x * (1 + g) * A * sigY
      dT = y * (1 + g) * sigY / (g * alpha * E)

      T1, T2, period = cycle_generator(dT, T0, Tdot)
      res = drivers.twobar_test(model, g * A, A, T1, T2, period, P, load_time, ncycles,
          verbose = False)

      points.append([x,y])
      conditions.append(res['classification'])

  return points, conditions, osx, osy, nx, ny

def plot_diagram(g, dx, dy, nx, ny, conditions):
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
  xpoints = [0.0, (g-1)/(g+1), 1]
  xpixels = [xmap(x) for x in xpoints]
  
  plt.xticks(xpixels, xpoints)

  ymap = lambda y: ny - 1 - (y / dy - 0.5)
  ypoints = [0.0, 1.0, 2.0]
  ypixels = [ymap(y) for y in ypoints]

  plt.yticks(ypixels, ypoints)

  plt.show()

if __name__ == "__main__":
  gamma = 3.0

  points, conditions, dx, dy, nx, ny = grid_classical(gamma)
  
  plot_diagram(gamma, dx, dy, nx, ny, conditions)
