#!/usr/bin/env python

import numpy as np

from neml import models, parse, uniaxial, elasticity, surfaces
from neml.drivers import newton, classify
import matplotlib.pyplot as plt
from matplotlib import rc

import multiprocessing

class ThreeBarProblem(object):
  """
    Two bars in series, one in parallel.

    This is the simplified model of the test specimen.
  """
  def __init__(self, l1, l2, l3, A1, A2, A3, mat1, mat2, mat3, 
      T0):
    """
      Parameters:
        l1      length of bar 1
        l2      length of bar 2
        l3      length of bar 3
        A1      area of bar 1
        A2      area of bar 2
        A3      area of bar 3
        mat1    material model for bar 1
        mat2    material model for bar 2
        mat3    material model for bar 3
        T0      initial temperature

      Notes:
        l1 = 2 * l2 + l3 or l3 = l1 - 2 * l2
    """
    self.l1 = l1
    self.l2 = l2
    self.l3 = l3

    self.A1 = A1
    self.A2 = A2
    self.A3 = A3

    self.mat1 = uniaxial.UniaxialModel(mat1)
    self.mat2 = uniaxial.UniaxialModel(mat2)
    self.mat3 = uniaxial.UniaxialModel(mat3)

    self.s1 = [0.0]
    self.s2 = [0.0]
    self.s3 = [0.0]

    self.e1 = [0.0]
    self.e2 = [0.0]
    self.e3 = [0.0]

    self.em1 = [0.0]
    self.em2 = [0.0]
    self.em3 = [0.0]

    self.et1 = [0.0]
    self.et2 = [0.0]
    self.et3 = [0.0]

    self.h1 = [self.mat1.init_store()]
    self.h2 = [self.mat2.init_store()]
    self.h3 = [self.mat3.init_store()]

    self.T = [T0]
    self.t = [0.0]

    self.u = [0.0]
    self.p = [0.0]

  def take_step(self, dt, T, verbose = False):
    """
      Step the system by incrementing time by dt and changing temperature to T

      Parameters:
        dt      time step
        T       next temperature
    """
    t_np1 = self.t[-1] + dt

    x0 = np.array([self.e1[-1], self.e2[-1], self.e3[-1]])

    def RJ(x):
      R = np.zeros((3,))
      J = np.zeros((3,3))
      
      a1 = self.mat1.alpha(T)
      a1p = self.mat1.alpha(self.T[-1])
      et1 = self.et1[-1] + (T - self.T[-1]) * (a1 + a1p) / 2.0 
      s1, h1, A1, u1, p1 = self.mat1.update(x[0] - et1, self.em1[-1], T, 
          self.T[-1], t_np1, self.t[-1], self.s1[-1], self.h1[-1], 0.0, 0.0)

      a2 = self.mat2.alpha(T)
      a2p = self.mat2.alpha(self.T[-1])
      et2 = self.et2[-1] + (T - self.T[-1]) * (a2 + a2p) / 2.0 
      s2, h2, A2, u2, p2 = self.mat2.update(x[1] - et2, self.em2[-1], T, 
          self.T[-1], t_np1, self.t[-1], self.s2[-1], self.h2[-1], 0.0, 0.0)

      a3 = self.mat3.alpha(T)
      a3p = self.mat3.alpha(self.T[-1])
      et3 = self.et3[-1] + (T - self.T[-1]) * (a3 + a3p) / 2.0 
      s3, h3, A3, u3, p3 = self.mat3.update(x[2] - et3, self.em3[-1], T, 
          self.T[-1], t_np1, self.t[-1], self.s3[-1], self.h3[-1], 0.0, 0.0)

      R[0] = s2 * self.A2 - s3 * self.A3
      R[1] = s2 * self.A2 + s1 * self.A1
      R[2] = x[0] * self.l1 - x[1] * self.l2 - x[2] * self.l3

      J[0,0] = 0.0
      J[0,1] = A2 * self.A2
      J[0,2] = -A3 * self.A3

      J[1,0] = A1 * self.A1
      J[1,1] = A2 * self.A2
      J[1,2] = 0.0

      J[2,0] = self.l1
      J[2,1] = -self.l2
      J[2,2] = -self.l3

      return R, J

    x = newton(RJ, x0, verbose = verbose)

    a1 = self.mat1.alpha(T)
    a1p = self.mat1.alpha(self.T[-1])
    et1 = self.et1[-1] + (T - self.T[-1]) * (a1 + a1p) / 2.0 
    s1, h1, A1, u1, p1 = self.mat1.update(x[0] - et1, self.em1[-1], T, 
        self.T[-1], t_np1, self.t[-1], self.s1[-1], self.h1[-1], 0.0, 0.0)

    a2 = self.mat2.alpha(T)
    a2p = self.mat2.alpha(self.T[-1])
    et2 = self.et2[-1] + (T - self.T[-1]) * (a2 + a2p) / 2.0 
    s2, h2, A2, u2, p2 = self.mat2.update(x[1] - et2, self.em2[-1], T, 
        self.T[-1], t_np1, self.t[-1], self.s2[-1], self.h2[-1], 0.0, 0.0)

    a3 = self.mat3.alpha(T)
    a3p = self.mat3.alpha(self.T[-1])
    et3 = self.et3[-1] + (T - self.T[-1]) * (a3 + a3p) / 2.0 
    s3, h3, A3, u3, p3 = self.mat3.update(x[2] - et3, self.em3[-1], T, 
        self.T[-1], t_np1, self.t[-1], self.s3[-1], self.h3[-1], 0.0, 0.0)

    self.s1.append(s1)
    self.s2.append(s2)
    self.s3.append(s3)

    self.e1.append(x[0])
    self.e2.append(x[1])
    self.e3.append(x[2])

    self.et1.append(et1)
    self.et2.append(et2)
    self.et3.append(et3)

    self.em1.append(x[0] - et1)
    self.em2.append(x[1] - et2)
    self.em3.append(x[2] - et3)

    self.h1.append(h1)
    self.h2.append(h2)
    self.h3.append(h3)

    self.t.append(self.t[-1] + dt)
    self.T.append(T)

    self.u.append(self.u[-1] + self.l1 * self.A1 * u1 +
        self.l2 * self.A2 * u2 + self.l3 * self.A3 * u3)
    self.p.append(self.p[-1] + self.l1 * self.A1 * p1 + 
        self.l2 * self.A2 * p2 + self.l3 * self.A3 * p3)

  def run_standard_cycle(self, T1, T2, t1, t2, th, nsteps1, nsteps2, 
      nstepsh, ncycles, verbose = False):
    """
      Ramp from T1 to T2 over t1, hold for th, ramp from T2 to T1 over t2,
      repeat ncycles times

      Parameters:
        T1          starting cycle temperature
        T2          ending cycle temperature
        t1          ramp time from T1 to T2
        t2          ramp time from T2 to T1
        th          hold time at T2
        nsteps1     number of steps in ramp up
        nsteps2     number of steps in ramp down
        nstepsh     number of steps during hold
        ncycles     number of times to repeat

      Optional:
        verbose     print a lot of info
    """
    dT1 = (T2 - T1) / nsteps1
    dt1 = t1 / nsteps1
    if nstepsh == 0:
      dth = 0
    else:
      dth = th / nstepsh
    dT2 = (T1 - T2) / nsteps2
    dt2 = t2 / nsteps2
    
    Tc = T1
    for i in range(ncycles):
      for j in range(nsteps1):
        Tc += dT1
        self.take_step(dt1, Tc, verbose = verbose)

      for j in range(nstepsh):
        self.take_step(dth, Tc, verbose = verbose)

      for j in range(nsteps2):
        Tc += dT2
        self.take_step(dt2, Tc, verbose = verbose)

def elastic_model(E, a, nu = 0.3):
  Emodel = elasticity.YoungsModulus(E)
  vmodel = elasticity.PoissonsRatio(nu)
  elmodel = elasticity.IsotropicLinearElasticModel(Emodel, vmodel)

  return models.SmallStrainElasticity(elmodel, alpha = a)

def verify():
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

  mat1 = elastic_model(E1, a1)
  mat2 = elastic_model(E2, a2)
  mat3 = elastic_model(E3, a3)

  print("Expected: %f" % expected)

  model = ThreeBarProblem(l1, l2, l3, A1, A2, A3, mat1, mat2, mat3, T1)
  model.run_standard_cycle(T1, T2, 1.0, 1.0, 0, 1, 1, 0, 1, verbose = False)
  calced = model.em3[1]

  print("Calculated: %f" % calced)

def compare_gr91_316_elastic():
  mat_inner = parse.parse_xml('gr91.xml', 'gr91_elastic')
  mat_outer = parse.parse_xml('316.xml', '316_elastic')
  
  # Remember we take the 1/2 symmetry view
  l1 = 58.7375
  l2 = 46.0375
  l3 = 12.7

  A1 = np.pi*(19.05**2.0 - (19.05-6.35)**2.0)
  A2 = np.pi*9.525**2.0 
  A3 = np.pi*3.175**2.0
 
  T1 = 160.0 + 273.15
  T2 = 548.9 + 273.15
  model_1 = ThreeBarProblem(l1, l2, l3, A1, A2, A3, mat_outer, mat_inner,
      mat_inner, T1)
  model_1.run_standard_cycle(T1, T2, 1.0, 1.0, 0, 10, 10, 0, 1, verbose = False)
  range_1 = np.max(model_1.em3) - np.min(model_1.em3)

  T1 = 354.4 + 273.15
  T2 = 548.9 + 273.15
  model_2 = ThreeBarProblem(l1, l2, l3, A1, A2, A3, mat_outer, mat_inner,
      mat_inner, T1)
  model_2.run_standard_cycle(T1, T2, 1.0, 1.0, 0, 10, 10, 0, 1, verbose = False)
  range_2 = np.max(model_2.em3) - np.min(model_2.em3)
  
  print("320 F to 1020 F: %f" % range_1)
  print("670 F to 1020 F: %f" % range_2)

def run_gr91_316_inelastic():
  mat_inner = parse.parse_xml('gr91.xml', 'gr91_simple')
  mat_outer = parse.parse_xml('316.xml', '316_simple')
  
  # Remember we take the 1/2 symmetry view
  l1 = 58.7375
  l2 = 46.0375
  l3 = 12.7

  A1 = np.pi*(19.05**2.0 - (19.05-6.35)**2.0)
  A2 = np.pi*9.525**2.0 
  A3 = np.pi*3.175**2.0
  
  """
  T1 = 160.0 + 273.15
  T2 = 548.9 + 273.15
  ramp = 1*3600.0
  hold = 1000*3600.0
  model_1 = ThreeBarProblem(l1, l2, l3, A1, A2, A3, mat_outer, mat_inner,
      mat_inner, T1)
  model_1.run_standard_cycle(T1, T2, ramp, ramp, hold, 25, 25, 25, 2, 
      verbose = True)
  """

  T1 = 354.4 + 273.15
  T2 = 548.9 + 273.15
  ramp = 1*3600.0
  hold = 1000*3600.0
  model_2 = ThreeBarProblem(l1, l2, l3, A1, A2, A3, mat_outer, mat_inner,
      mat_inner, T1)
  model_2.run_standard_cycle(T1, T2, ramp, ramp, hold, 25, 25, 25, 20, 
      verbose = False)
  
  print("Strain range: %f" % (np.max(model_2.em3[-75:]) - 
    np.min(model_2.em3[-75:])))

  plt.plot(np.array(model_2.t)/3600.0, model_2.em3)
  plt.xlabel("Time (hrs)")
  plt.ylabel("Mechanical strain in gauge bar (mm/mm)")
  plt.show()

def plastic_model(E, Y, a, nu = 0.3):
  Emodel = elasticity.YoungsModulus(E)
  vmodel = elasticity.PoissonsRatio(nu)
  elmodel = elasticity.IsotropicLinearElasticModel(Emodel, vmodel)
  surf = surfaces.IsoJ2()
  
  return models.SmallStrainPerfectPlasticity(elmodel, surf, Y, alpha = a)

def bree_point(A1, A2, A3, l1, l2, l3, E1, E2, E3, Y1, Y2, Y3, a1, a2, a3, 
    T1, T2, ncycles = 5, nsteps = 25):
  mat1 = plastic_model(E1, Y1, a1)
  mat2 = plastic_model(E2, Y2, a2)
  mat3 = plastic_model(E3, Y3, a3)

  model = ThreeBarProblem(l1, l2, l3, A1, A2, A3, mat1, mat2, mat3, T1)
  nsteps_total = 2 * nsteps

  try:
    model.run_standard_cycle(T1, T2, 1.0, 1.0, 0, nsteps, nsteps, 0,
        ncycles, verbose = False)
  except Exception:
    return "collapse"

  ua = model.u[-nsteps_total-1]
  ub = model.u[-1]

  pa = model.p[-nsteps_total-1]
  pb = model.p[-1]

  e1a = model.e1[-nsteps_total-1]
  e1b = model.e1[-1]

  e2a = model.e3[-nsteps_total-1]
  e2b = model.e3[-1]
  
  return classify(ua, ub, pa, pb, e1a, e1b, e2a, e2b)

def make_grid(xmax, ymax, nx, ny):
  xp = np.linspace(0, xmax, nx+1)
  xc = np.diff(xp) / 2 + xp[:-1]
  yp = np.linspace(0, ymax, ny+1)
  yc = np.diff(yp) / 2 + yp[:-1]
  
  points = []
  for x in xc:
    for y in yc:
      points.append([x,y])

  return np.array(points)

def plot_diagram(x, y, nx, ny, conditions, xlabel = None, ylabel = None):
  dx = x / nx
  dy = y / ny

  plt.figure()
  rc('text', usetex = True)
  rc('font', size = 18)

  cmap = {
      "collapse": (51.0/255,160/255.0,44.0/255),
      "elastic" : (166.0/255,206.0/255,227.0/255),
      "elastic shakedown": (31.0/255,120.0/255,180.0/255),
      "plastic shakedown": (178.0/255,223.0/255,138.0/255),
      "ratcheting": (51.0/255,160/255.0,44.0/255)
      }
  img = np.zeros((ny, nx, 3))
  k = 0
  for i in range(nx):
    for j in range(ny):
      img[j,i] = cmap[conditions[k]]
      k += 1

  img = img[::-1,:,:]

  plt.imshow(img, interpolation = 'none', aspect = 'equal')
  
  xmap = lambda x: x / dx - 0.5
  xpoints = np.linspace(0,x,5)

  xpixels = [xmap(x) for x in xpoints]
  plt.xticks(xpixels, xpoints)

  ymap = lambda y: ny - 1 - (y / dy - 0.5)
  ypoints = np.linspace(0,y,5)
  ypixels = [ymap(y) for y in ypoints]
  plt.yticks(ypixels, ypoints)

  plt.xlim([xmap(0.0), xmap(x)])
  if xlabel is not None:
    plt.xlabel(xlabel)
  if ylabel is not None:
    plt.ylabel(ylabel)
  plt.tight_layout()
  plt.show()

if __name__ == "__main__":
  verify()
  compare_gr91_316_elastic()
  run_gr91_316_inelastic()

  l1 = 58.7375
  l2 = 46.0375
  l3 = 12.7

  A1 = np.pi*(19.05**2.0 - (19.05-6.35)**2.0)
  A2 = np.pi*9.525**2.0 
  A3 = np.pi*3.175**2.0

  T1 = 0

  Einner = 190680.0
  ainner = 1.3016e-5
  Yinner = 370.2

  Eouter = 171760.0
  aouter = 1.9316e-5/100
  Youter = 126.7
  
  def eval_l3l1(pt):
    x = pt[0]
    y = pt[1]
    l3 = x * l1
    T2 = y * Yinner / (ainner * Einner)
    return bree_point(A1, A2, A3, l1, l1 - l3, l3, Eouter, Einner, Einner, 
        Youter, Yinner, Yinner, aouter, ainner, ainner, T1, T2)
  
  x = 1.0
  y = 4.0
  nx = 50
  ny = 50
  grid = make_grid(x, y, nx, ny)

  pool = multiprocessing.Pool()
  conditions = pool.map(eval_l3l1, grid)
  plot_diagram(x, y, nx, ny, conditions, xlabel = r'$\frac{l_3}{l_1}$', 
      ylabel = r'$\frac{\alpha_{inner} E_{inner}\Delta T}{Y_{inner}}$')
  
  def eval_A3A1(pt):
    x = pt[0]
    y = pt[1]
    A3 = x * A2
    T2 = y * Yinner / (ainner * Einner)
    return bree_point(A1, A2, A3, l1, l1 - l3, l3, Eouter, Einner, Einner, 
        Youter, Yinner, Yinner, aouter, ainner, ainner, T1, T2)
  
  x = 2.0
  y = 4.0
  nx = 50
  ny = 50
  grid = make_grid(x, y, nx, ny)

  pool = multiprocessing.Pool()
  conditions = pool.map(eval_A3A1, grid)
  plot_diagram(x, y, nx, ny, conditions, xlabel = r'$\frac{A_3}{A_1}$', 
      ylabel = r'$\frac{\alpha_{inner} E_{inner}\Delta T}{Y_{inner}}$')

  def eval_A2A1(pt):
    x = pt[0]
    y = pt[1]
    A2 = x * A1
    T2 = y * Yinner / (ainner * Einner)
    return bree_point(A1, A2, A3, l1, l1 - l3, l3, Eouter, Einner, Einner, 
        Youter, Yinner, Yinner, aouter, ainner, ainner, T1, T2)

  x = 2.0
  y = 4.0
  nx = 50
  ny = 50
  grid = make_grid(x, y, nx, ny)

  pool = multiprocessing.Pool()
  conditions = pool.map(eval_A2A1, grid)
  plot_diagram(x, y, nx, ny, conditions, xlabel = r'$\frac{A_2}{A_1}$', 
      ylabel = r'$\frac{\alpha_{inner} E_{inner}\Delta T}{Y_{inner}}$')
