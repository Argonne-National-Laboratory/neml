#!/usr/bin/env python

import numpy as np
import numpy.linalg as la
from numpy.polynomial.legendre import leggauss as lgg
import scipy.optimize as opt

def generate_thickness_gradient(ri, ro, T1, T2, Tdot_hot, hold, 
    Tdot_cold = None, hold_together = 0.0, delay = 0.0):
  """
    Generate our usual through-wall temperature gradient

    Parameters:
      ri                inner radius
      ro                outer radius
      T1                initial temperature
      T2                hot temperature
      Tdot_hot          heating rate, also cooling rate if not specified
      hold              hold at gradient

    Optional:
      Tdot_cold         cooling rate
      hold_together     hold at no temperature gradient
  """
  if Tdot_cold is None:
    Tdot_cold = Tdot

  dT = T2 - T1
  thot = np.abs(dT / Tdot_hot)
  tcold = np.abs(dT / Tdot_cold)
  period = thot + hold + tcold + hold_together

  def gradient(r, t):
    if t < delay:
      return T1
    t = (t-delay) % period

    if t < thot:
      Tright = T1 + t * Tdot_hot * np.sign(dT)
    elif t < thot + hold:
      Tright = T2
    elif t < thot + hold + tcold:
      Tright = T2 - (t - thot - hold) * Tdot_cold * np.sign(dT)
    else:
      return T1

    xi = (r - ri) / (ro - ri)

    return (1-xi) * T1 + xi * Tright

  return gradient

def mesh(rs, ns, bias = False, n = 2.0):
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
      bias = False, factor = 2.0, ngpts = 2):
    """
      Parameters:
        rs      radii deliminating each region
        mats    material model for each region
        ns      number of elements per region
        T       temperature as a function of (r, t)
        factor  pressure, as a function of (t)

      Optional:
        rtol    relative tolerance for N-R solve
        atol    absolute tolerance for N-R solve
        bias    if true bias the mesh to the edges
        n       bias factor
        ngpts   number of gauss points to use per element
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

    self.mesh = mesh(self.rs, self.ns, bias = bias, n = factor)
    self.nnodes = len(self.mesh)
    self.nelem = len(self.mesh) - 1

    self.materials = [m for m,n in zip(self.mats, self.ns) for i in range(n)]
    
    self.times = [0.0]
    self.pressures = [self.p(0)]
    self.temperatures = [np.array([self.T(xi, 0) for xi in self.mesh])]
    
    self.displacements = [[0.0 for i in range(self.nnodes)]]
    self.axialstrain = [0.0]
    self.stresses = [[np.zeros((4,)) for n in self.ns for i in range(n)]]
    self.strains = [[np.zeros((4,)) for n in self.ns for i in range(n)]]
    self.mstrains = [[np.zeros((4,)) for n in self.ns for i in range(n)]]
    self.estrains = [[np.zeros((4,)) for n in self.ns for i in range(n)]]

    self.gpoints, self.gweights = lgg(ngpts)

    self.histories = [[[m.init_store() for xi in self.gpoints]
      for m in self.materials]]

    self.ls = np.diff(self.mesh)

    # Save a bit of time
    self.ri = np.array([[(xi+1)/2*li + ri for xi in self.gpoints] for li, ri 
      in zip(self.ls, self.mesh)])
    self.Nl = np.array([[(1-xi)/2, (1+xi)/2] for xi in self.gpoints])
    self.Bl = np.array([[-1, 1] for xi in self.gpoints])

  def step(self, t):
    """
      Update to the next state

      Parameters:
        t       time requested
    """
    T = np.array([self.T(xi, t) for xi in self.mesh])
    p = self.p(t)
    dt = t - self.times[-1]
    
    # Get a guess
    x0 = np.concatenate((self.displacements[-1], [self.axialstrain[-1]]))

    R0 = self.R(x0, T, p)

    self.times.append(t)
    self.pressures.append(p)
    self.temperatures.append(T)

  def R(self, x, T, p):
    """
      Compute the residual and the updated history

      Parameters:
        x       iterate
        T       nodal temperatures
        p       left pressure
    """
    d = x[:-1]
    ez = x[-1]
    
    Fext = np.zeros((self.nnodes+1,))
    Fext[0] = self.rs[-1] * p
    Fext[-1] = p * self.r / self.t

    Fint = np.zeros((self.nnodes+1,))
    for e, (mat, l, rs) in enumerate(zip(self.materials, self.ls, self.ri)):
      de = d[e:e+2]
      Te = T[e:e+2]
      Tp = self.temperatures[-1][e:e+2]
      dT = Te - Tp
      err = np.dot(self.Bl, de) * 2.0 / l
      ett = np.dot(self.Nl, de) / rs
      ezz = np.array([ez, ez])
      strain = np.vstack((err, ett, ezz)).T

    return Fint - Fext
