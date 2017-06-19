#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, neml, elasticity, drivers, surfaces, hardening, ri_flow, uniaxial

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

import multiprocessing

def single_material_generator(E, nu, sY, alpha, n):
  models = []
  for i in range(n):
    youngs = elasticity.YoungsModulus(E)
    poisson = elasticity.PoissonsRatio(nu)
    elastic = elasticity.IsotropicLinearElasticModel(youngs, poisson)
    surface = surfaces.IsoJ2()
    models.append(
        uniaxial.UniaxialModel(neml.SmallStrainPerfectPlasticity(
          elastic, surface, sY, alpha = alpha)))

  return models

def classify_case(strain, energy, work, ncycle, rtol = 1.0e-4, atol = 1.0e-10):
  ub = energy[-1]
  ua = energy[-1-ncycle]
  pb = work[-1]
  pa = work[-1-ncycle]
  eb = strain[-1]
  ea = strain[-1-ncycle]
  
  if np.abs(pb) < atol:
    return 'elastic'
  elif np.abs(ub-ua) < rtol * np.abs(ub):
    return 'elastic shakedown'
  elif np.abs(eb-ea) < rtol * np.abs(eb):
    return 'plastic shakedown'
  else:
    return 'ratcheting'

def run_case(P, dT, models, nsteps_cycle = 25):
  try:
    time, strain, energy, work = drivers.bree(models,  P, dT,
        nsteps_cycle = 25)
  except Exception:
    return "collapse"
  return classify_case(strain, energy, work, nsteps_cycle * 2)

def plot_diagram(dx, dy, nx, ny, conditions):
  """
    Makes a pretty plot of the results
  """
  rc('text', usetex = True)

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

  plt.imshow(img, interpolation = 'none', aspect = 0.3)
  
  xmap = lambda x: x / dx - 0.5
  xpoints = [0.0, 0.5, 1]
  xpixels = [xmap(x) for x in xpoints]
  
  plt.xticks(xpixels, xpoints)

  ymap = lambda y: ny - 1 - (y / dy - 0.5)
  ypoints = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  ypixels = [ymap(y) for y in ypoints]

  plt.yticks(ypixels, ypoints)

  plt.plot([xmap(0.0), xmap(1.0)], [ymap(1.0), ymap(0.0)], 'k-')
  plt.plot([xmap(0.0), xmap(0.5)], [ymap(2.0), ymap(2.0)], 'k-')
  plt.plot([xmap(0.5), xmap(1.0)], [ymap(2.0), ymap(0.0)], 'k-')

  yvs = np.linspace(2.0, 5.0)
  xvs = 1.0 / yvs

  plt.plot([xmap(x) for x in xvs], [ymap(y) for y in yvs], 'k-')

  plt.xlim([xmap(0.0), xmap(1.0)])
  plt.xlabel(r'$\frac{p}{\sigma_y}$')
  plt.ylabel(r'$\frac{\alpha E \Delta T}{\sigma_y}$')

  plt.show()

def runme(p):
  x = p[0]
  y = p[1]
  print("(%f,%f)" % (x,y))
  P = sY * x
  dT = sY * y * 2 / (E * alpha)
  return run_case(P, dT, models)

def grid_classical(models, E, nu, sY, alpha, 
    max_x = 1.1, max_y = 4.0, nx = 11, ny = 20):
  osx = max_x / (nx)
  osy = max_y / (ny)
  
  points = []
  for x in np.linspace(osx/2, max_x-osx/2, nx):
    for y in np.linspace(osy/2, max_y-osy/2, ny):
      points.append([x,y])

  pool = multiprocessing.Pool()
  conditions = pool.map(runme, points)

  return points, conditions, osx, osy, nx, ny

if __name__ == "__main__":
  E = 100000.0
  nu = 0.25
  sY = 100.0
  alpha = 1.0e-5
  n = 4

  models = single_material_generator(E, nu, sY, alpha, n)

  points, conditions, dx, dy, nx, ny = grid_classical(models,
      E, nu, sY, alpha, max_x = 0.9999, max_y = 5.0, nx = 15, ny = 75)

  plot_diagram(dx, dy, nx, ny, conditions)

