#!/usr/bin/env python3

import sys
sys.path.append('../../..')

from neml.math import rotations
from neml.cp import harmonics, polefigures, crystallography

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def sphere_dist(x1,x2):
  d = np.dot(x1,x2)
  if d > 1.0:
    d = 1.0
  elif d < -1.0:
    d = -1.0
  return np.arccos(d)

if __name__ == "__main__":
  N = 100

  pts, wts = harmonics.quadrature_S1(N)
  
  plt.figure()
  plt.scatter(np.cos(pts), np.sin(pts))
  plt.axes().set_aspect('equal', 'datalim')
  plt.show()
  
  dds = []
  for n in range(1,N):
    pts, wts = harmonics.quadrature_S1(n)
    nn = len(pts)
    dist = np.zeros((nn-1,))
    for i,pt in enumerate(pts[1:]):
      dist[i] = np.abs(pts[0]-pts[i+1]) % (2.0 * np.pi)
    dds.append(np.min(dist))

  plt.loglog(np.array(range(1,N))+1,dds)
  plt.show()

  N = 10

  theta, phi, wts = harmonics.quadrature_S2(N)

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta),
      zs = np.cos(phi))
  plt.xlim([-1,1])
  plt.ylim([-1,1])
  ax.set_aspect('equal')
  plt.show()

  N = 50
  
  npts = []
  dds = []
  for n in range(1,N):
    theta, phi, wts = harmonics.quadrature_S2(n)
    npts.append(len(phi))
    pts = np.array([np.sin(phi) * np.cos(theta), np.sin(phi)*np.sin(theta), 
      np.cos(phi)]).T
    pt0 = pts[0]
    dds.append(np.inf)
    for pt in pts[1:]:
      di = sphere_dist(pt0,pt)
      if di < dds[-1]:
        dds[-1] = di
  
  plt.loglog(npts, dds)
  plt.show()

  N = 5
  pts, wts = harmonics.quadrature_SO3(N)
  lattice = crystallography.CubicLattice(1.0)
  lattice.add_slip_system([1,1,0],[1,1,1])
  polefigures.pole_figure_discrete(pts, [1,1,1],lattice)
  plt.show()

  N = 15
  npts = []
  dds = []
  for n in range(1,N):
    pts, wts = harmonics.quadrature_SO3(n)
    npts.append(len(pts))
    dds.append(np.inf)
    for pt in pts[1:]:
      d = pts[0].distance(pt)
      if d < dds[-1]:
        dds[-1] = d

  plt.loglog(npts, dds)
  plt.show()
