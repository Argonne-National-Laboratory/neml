#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import axisym, models, elasticity, surfaces, creep

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import scipy.interpolate as inter

def progress(f, l = 20):
  """
    Progress bar

    Parameters:
      f         completed fraction

    Optional:
      l         length
  """
  fill = int(round(f*l))
  pprint = "\rProgress: [%s] %6.2f%%" % ("#"*fill + "-"*(l-fill), f*100)
  sys.stdout.write(pprint)
  sys.stdout.flush()

def gen_material(E, nu, sY, alpha, A, n):
  """
    Generate a perfectly plastic + creep material model.
  """
  elastic = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, 
      "poissons")
  surface = surfaces.IsoJ2()

  pmodel = models.SmallStrainPerfectPlasticity(elastic, surface, sY)
  smodel = creep.PowerLawCreep(A, n)
  cmodel = creep.J2CreepModel(smodel)

  return models.SmallStrainCreepPlasticity(elastic, pmodel, cmodel, alpha = alpha)

if __name__ == "__main__":
  base = gen_material(156000.0, 0.31, 102, 1.99e-5, 6.94e-22, 6.0)
  clad = gen_material(174000.0, 0.3, 43.0, 1.58e-5, 4.66e-19, 4.6)

  Ro = 1750.0
  t = 60.0
  tclad = 0.6
  
  p = 1.74
  tramp = 5.0 * 3600.0
  nramp = 25

  theat = 5.0 * 3600.0
  nheat = 50
  thold = 10000.0 * 3600.0
  nhold = 5
  tcool = 5.0 * 3600.0
  ncool = 50

  ncycles = 30

  msize = 0.06
  bsize = 1.0

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
      [clad, base], [int(tclad / msize),int((t-tclad)/bsize)],
      gradient, pressure, bias = True)

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
    amodel.step(t)
    progress(float(i)/(len(times)-1))
  print("")
  
  # These are working back from 
  stresses = np.array(amodel.stresses)
  mstrains = np.array(amodel.mstrains)
  gpoints = np.array(amodel.ri)

  # Get element averages?
  epoints = np.mean(gpoints, axis = 1)
  estresses = np.mean(stresses, axis = 2)
  emstrains = np.mean(mstrains, axis = 2)
   

