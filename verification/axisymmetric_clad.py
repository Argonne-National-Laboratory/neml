#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import axisym, neml, elasticity, surfaces, creep

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import scipy.interpolate as inter

def gen_material(E, nu, sY, alpha, A, n):
  """
    Generate a perfectly plastic + creep material model.
  """
  youngs = elasticity.YoungsModulus(E)
  poissons = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(youngs, poissons)
  surface = surfaces.IsoJ2()

  pmodel = neml.SmallStrainPerfectPlasticity(elastic, surface, sY)
  smodel = creep.PowerLawCreep(A, n)
  cmodel = creep.J2CreepModel(smodel)

  return neml.SmallStrainCreepPlasticity(pmodel, cmodel, alpha = alpha)

if __name__ == "__main__":
  vdata = np.loadtxt('warp-3d-data.csv', delimiter = ',', skiprows = 1)

  base = gen_material(156000.0, 0.31, 102, 1.99e-5, 6.94e-22, 6.0)
  clad = gen_material(174000.0, 0.3, 43.0, 1.58e-5, 4.66e-19, 4.6)

  Ro = 1750.0
  t = 60.0
  tclad = 3.0
  
  p = 0.85
  tramp = 5.0 * 3600.0
  nramp = 25

  theat = 5.0 * 3600.0
  nheat = 50
  thold = 10000.0 * 3600.0
  nhold = 5
  tcool = 5.0 * 3600.0
  ncool = 50

  ncycles = 30

  msize = 1.0

  Ri = Ro - t

  T1 = 0.0
  T2 = 90.0

  Tdot_hot = (T2 - T1) / theat
  Tdot_cold = (T2 - T1) /tcool
  ntotal = nheat + nhold + ncool

  gradient = axisym.generate_thickness_gradient(Ro, Ri, T1, T2,
      Tdot_hot, thold, Tdot_cold = Tdot_cold, delay = tramp)

  def pressure(t):
    if t < tramp:
      return p/tramp * t
    else:
      return p

  amodel = axisym.AxisymmetricProblem([Ri, Ri+tclad, Ro], 
      [clad, base], [int(tclad / msize),int((t-tclad)/msize)],
      gradient, pressure)

  ts_period = np.concatenate(
      (
        np.linspace(0, theat, num = nheat),
        np.linspace(theat, theat + thold, num = nhold),
        np.linspace(theat + thold, theat + thold + tcool, num = ncool)
        ))

  ts_ramp = np.linspace(0, tramp, nramp)

  times = [ts_ramp]
  for i in range(ncycles):
    times.append(ts_period + times[-1][-1])
  
  times = np.concatenate(times)

  for i,t in enumerate(times):
    print(i)
    amodel.step(t, verbose = True)

  # These are working back from 
  stresses = np.array(amodel.stresses)
  mstrains = np.array(amodel.mstrains)
  gpoints = np.array(amodel.ri)

  # Get element averages?
  epoints = np.mean(gpoints, axis = 1)
  estresses = np.mean(stresses, axis = 2)
  emstrains = np.mean(mstrains, axis = 2)
   
  # Order rr, tt, zz, rz
  cols_me = [0,1,2,5]
  cols_warp = [28,29,30,33]
  cols_moose = [6,7,8,9]
  colors = ['k', 'r', 'b', 'g']
  for cm,cw,cmo,c in zip(cols_me, cols_warp, cols_moose, colors):
    plt.plot(epoints, emstrains[-1,:,cm], color = c, ls = '-')
    plt.plot(vdata[:,65], vdata[:,cw], color = c, ls = '--')

  plt.xlabel("Position (mm)")
  plt.ylabel("Strain (mm/mm)")
  plt.tight_layout()

  plt.show()

  cols_warp = [2,3,4,7]
  cols_moose = [10,11,12,13]
  for cm,cw,cmo,c in zip(cols_me, cols_warp, cols_moose, colors):
    plt.plot(epoints, estresses[-1,:,cm], color = c, ls = '-')
    plt.plot(vdata[:,65], vdata[:,cw], color = c, ls = '--')
  
  plt.xlabel("Position (mm)")
  plt.ylabel("Stress (MPa)")
  plt.tight_layout()
  plt.show()
