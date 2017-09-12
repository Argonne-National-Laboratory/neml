#!/usr/bin/env python

import numpy as np
import numpy.linalg as la
import scipy.optimize as opt
from drivers import newton

import matplotlib.pyplot as plt

def mesh(rs, ns, bias = False, n = 4):
  if not bias:
    xpoints = [np.linspace(rs[0], rs[1], ns[0]+1)]
    for i in range(1, len(ns)):
      xpoints.append(np.linspace(rs[i], rs[i+1], ns[i]+1)[1:])
    xpoints = np.concatenate(tuple(xpoints))
  else:
    xpoints = [rs[0]]
    for i in range(len(rs)-1):
      if ns[i] % 2 == 0:
        even = True
        nuse = ns[i] / 2 + 1
      else:
        even = False
        nuse = (ns[i] + 1) / 2 + 1

      x1 = (np.linspace(0,1, nuse)[1:])**n
      x2 = 1-(np.linspace(0,1, nuse)[:-1])**n + 1.0
      if even:
        xss = np.concatenate((x1, np.flip(x2,0)))/2
      else:
        xss = np.concatenate((x1[:-1], np.flip(x2,0)))/2
      
      # Transform
      xss = xss * (rs[i+1] - rs[i]) + rs[i]
      xpoints.extend(list(xss))
    
    xpoints = np.array(xpoints)

  return xpoints

class AxisymmetricProblem(object):
  """
    Python driver for our axisymmetric reduced problem.
  """
  def __init__(self, rs, mats, ns, T, p, rtol = 1.0e-6, atol = 1.0e-8,
      bias = False, n = 2.0):
    """
      Parameters:
        rs      radii deliminating each region
        mats    material model for each region
        ns      number of elements per region
        T       temperature as a function of (r, t)
        p       pressure, as a function of (t)

      Optional:
        rtol    relative tolerance for N-R solve
        atol    absolute tolerance for N-R solve
        bias    if true bias the mesh to the edges
        n       bias factor
    """
    # Check
    if len(rs) - 1 != len(mats):
      raise ValueError("Length of radii must be one more then length of regions")
    if len(mats) != len(ns):
      raise ValueError("Inconsistent region definitions")

    self.rs = rs
    self.mats = mats
    self.ns = ns
    self.T = T
    self.p = p

    self.rtol = rtol
    self.atol = atol

    self.ts = np.diff(rs)
    self.t = np.sum(self.ts)
    self.r = rs[-1]

    self.mesh = mesh(self.rs, self.ns, bias = bias, n = n)
    
    self.histories = []
    self.histories.append([[m.init_hist() for i in range(n) for m,n in zip(
      self.mats, self.ns)]])

    self.times = [0.0]
    self.pressures = [p(0)]
    self.temperatures = [[T(xi, 0) for xi in self.mesh]]

  def step(self, t):
    pass
