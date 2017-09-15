#!/usr/bin/env python

import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import scipy.linalg.lapack as lapack
from numpy.polynomial.legendre import leggauss as lgg
import scipy.optimize as opt

import nemlmath

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
    
    self.ngpts = ngpts
    self.gpoints, self.gweights = lgg(ngpts)

    self.histories = [[np.array([m.init_store() for xi in self.gpoints])
      for m in self.materials]]
    self.stresses = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]
    self.strains = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]
    self.tstrains = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]
    self.mstrains = [np.array([np.zeros((self.ngpts,6)) for n in self.ns for i in range(n)])]

    self.ls = np.diff(self.mesh)

    # Save a bit of time
    self.Nl = np.array([[(1-xi)/2, (1+xi)/2] for xi in self.gpoints])
    self.ri = np.array([np.dot(self.Nl, self.mesh[e:e+2]) for e in range(self.nelem)])
    self.Bl = np.array([[-1.0/2, 1.0/2] for xi in self.gpoints])

  def step(self, t, rtol = 1.0e-6, atol = 1.0e-8, verbose = False):
    """
      Update to the next state

      Parameters:
        t       time requested
    """
    T = np.array([self.T(xi, t) for xi in self.mesh])
    p = self.p(t)

    # Get a guess
    x = np.concatenate((self.displacements[-1], [self.axialstrain[-1]]))

    R, (J11_du,J11_d,J11_dl,J12,J21,J22), strains, tstrains, mstrains, stresses, histories = self.RJ(x, T, p, t)
    nR0 = la.norm(R)
    nR = nR0

    if verbose:
      print("Iter.\tNorm res")
      print("%i\t%e" % (0, nR0))
    
    i = 0
    while (nR > atol) and (nR / nR0 > rtol):
      R1 = R[:-1]
      R2 = R[-1]
      c = R1 - R2/J22 * J12
      u = -J12 / J22
      v = J21
      x1 = cute_solve(J11_du, J11_d, J11_dl, u, v, c)
      x2 = (R2 - np.dot(J21,x1)) / J22

      x[:-1] -= x1
      x[-1] -= x2
      
      R, (J11_du,J11_d,J11_dl,J12,J21,J22), strains, tstrains, mstrains, stresses, histories = self.RJ(x, T, p, t)
      nR = la.norm(R)
      i += 1

      if verbose:
        print("%i\t%e" % (i, nR))

    self.times.append(t)
    self.pressures.append(p)
    self.temperatures.append(T)
    self.strains.append(strains)
    self.tstrains.append(tstrains)
    self.mstrains.append(mstrains)
    self.stresses.append(stresses)
    self.histories.append(histories)

  def strain(self, d, l, r, ez):
    """
      Compute the full strains

      Parameters:
        d   displacements
        l   length
        r   radius
    """
    strain = np.zeros((len(r),6))
    strain[:,0] = np.dot(self.Bl, d) * 2.0 / l
    strain[:,1] = np.dot(self.Nl, d) / r
    strain[:,2] = ez 
    return strain

  def RJ(self, x, T, p, t):
    """
      Compute the residual, jacobian, and the rest of the updated quantities

      Parameters:
        x       iterate
        T       nodal temperatures
        Tn      previous temperatures
        p       left pressure
    """
    d = x[:-1]
    ez = x[-1]
    
    Fext = np.zeros((self.nnodes+1,))
    Fext[0] = p
    Fext[-1] = p * self.r / (2*self.t)
    
    strains = np.zeros((self.nelem, self.ngpts, 6))
    tstrains = np.zeros((self.nelem, self.ngpts, 6))
    mstrains = np.zeros((self.nelem, self.ngpts, 6))
    stresses = np.zeros((self.nelem, self.ngpts, 6))
    histories = []

    Fint = np.zeros((self.nnodes+1,))
    A11_dl = np.zeros((self.nnodes-1,))
    A11_d = np.zeros((self.nnodes,))
    A11_du = np.zeros((self.nnodes-1,))
    A12 = np.zeros((self.nnodes,))
    A21 = np.zeros((self.nnodes,))
    A22 = 0.0

    for e, (mat, l, rs) in enumerate(zip(self.materials, self.ls, self.ri)):
      de = d[e:e+2]
      T_e = np.dot(self.Nl, T[e:e+2])
      T_n = np.dot(self.Nl, self.temperatures[-1][e:e+2])
      
      therm_inc = np.array([mat.alpha(Ti)*(Ti-Tii)*np.array([1,1,1,0,0,0]) 
        for Ti, Tii in zip(T_e, T_n)])

      strain = self.strain(de, l, rs, ez)
      tstrain = self.tstrains[-1][e] + therm_inc
      mstrain = strain - tstrain
      mstrain_n = self.mstrains[-1][e]
      stress_n = self.stresses[-1][e]
      hist_n = self.histories[-1][e]

      stress = np.zeros((self.ngpts,6))
      history = np.zeros((self.ngpts,len(hist_n[0])))

      for i in range(self.ngpts):
        si, hi, ti, ui, pi = mat.update_sd(mstrain[i], mstrain_n[i],
            T_e[i], T_n[i], t, self.times[-1], stress_n[i], hist_n[i], 
            0.0, 0.0)
        stress[i] = si
        history[i] = hi
        
        wi = self.gweights[i]
        #Actually add in our terms
        Fint[e:e+2] += wi * (si[0] * self.Bl[i] * 2.0 / l + (si[1] - si[0]) * self.Nl[i] / rs[i]) * l / 2.0
        Fint[-1] += wi * si[2] * l / 2.0 / self.t
        
        DE = np.array([self.Bl[i] * 2.0/l, self.Nl[i] / rs[i]])
        
        # Remember wi and the jacobian
        J11 = np.dot(np.outer(2.0/l * self.Bl[i], ti[0,:2]) + np.outer(self.Nl[i] / rs[i], 
            ti[1,:2] - ti[0,:2]), DE) * wi * l / 2.0 
        J12 = (ti[0,2] * self.Bl[i] * 2.0 / l + (ti[1,2] - ti[0,2]) * self.Nl[i] / rs[i]) * wi * l / 2.0
        J21 = (np.dot(ti[2,:2], DE) / self.t) * wi * l / 2.0 
        J22 = (ti[2,2] / self.t) * wi * l / 2.0

        A11_du[e] += J11[0,1]
        A11_d[e] += J11[0,0]
        A11_d[e+1] += J11[1,1]
        A11_dl[e] += J11[1,0]
        A12[e:e+2] += J12
        A21[e:e+2] += J21
        A22 += J22

      strains[e] = strain
      tstrains[e] = tstrain
      mstrains[e] = mstrain
      stresses[e] = stress
      histories.append(history)

    return Fint - Fext, (A11_du,A11_d,A11_dl,A12,A21,A22), strains, tstrains, mstrains, stresses, histories


def cute_solve(DU, D, DL, u, v, c):
  """
    Do the cute solve of:

    (B + u v.T) x = c

    with B tri-diagonal
  """
  DL, D, DU, DU2, IPIV = nemlmath.dgttrf(DL, D, DU)
  
  x = nemlmath.dgttrs(DL,D,DU,DU2,IPIV, c)
  y = nemlmath.dgttrs(DL,D,DU,DU2,IPIV, u)
  
  return x - np.dot(v,x) / (1+np.dot(v,y)) * y
