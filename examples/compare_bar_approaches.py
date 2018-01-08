#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import neml, elasticity, surfaces, hardening, ri_flow, arbbar, drivers

import numpy as np
import matplotlib.pyplot as plt

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

def cycle_generator_delay(dT, T0, Tdot, delay):
  period = np.abs(2.0 * dT / Tdot)
  Tdot_sign = np.sign(dT) * Tdot

  T1 = lambda t: np.piecewise(t, 
      [t <= delay, 
        np.logical_and(t>delay, (t - delay) % period < period / 2),
        np.logical_and(t>delay, (t-delay) % period >= period / 2)],
      [lambda tt: (tt*0 + 1) * T0,
        lambda tt: T0 + Tdot_sign * ((tt-delay) % period),
        lambda tt: T0 + dT - Tdot_sign * (((tt-delay) % period) - period / 2)])

  T2 = lambda tt: T0 * (tt*0 + 1) 

  return T1, T2, period

def make_mat(E, nu, sy, k, a):
  youngs = elasticity.YoungsModulus(E)
  poissons = elasticity.PoissonsRatio(nu)
  elastic = elasticity.IsotropicLinearElasticModel(youngs, poissons)
  surface = surfaces.IsoJ2()
  harden = hardening.LinearIsotropicHardeningRule(sy, k)
  
  flow = ri_flow.RateIndependentAssociativeFlow(surface, harden)
  return neml.SmallStrainRateIndependentPlasticity(elastic, flow, alpha = a)

def old_way(A1, A2, mat, l, P, dT, ncycles, nsteps_load = 20, nsteps_cycle = 40):
  T1, T2, per = cycle_generator(dT, 0, 1)
  driver = drivers.Driver_sd_twobar(A1, A2, T1, T2, per, P, 1.0, mat,
      nsteps_load = nsteps_load, nsteps_cycle = nsteps_cycle)
  driver.load_up()
  for i in range(ncycles):
    driver.one_cycle()

  return driver.strain1_int[-1][0] * l

def new_way(A1, A2, mat, l, P, dT, ncycles, nsteps_load = 20, nsteps_cycle = 40):
  delay = 1
  T1, T2, per = cycle_generator_delay(dT, 0, 1, delay)

  ffn = lambda t: np.piecewise(t, [t < delay, t >= delay],
      [lambda tt: P / delay * tt, lambda tt: P])
  fixed = lambda t: 0.0
  
  model = arbbar.BarModel()
  
  model.add_node(1)
  model.add_node(2)

  model.add_edge(1,2,object=arbbar.Bar(mat, A1, l, T = T1))
  model.add_edge(1,2,object=arbbar.Bar(mat, A2, l, T = T2))

  times = np.linspace(0, delay, nsteps_load+1)[1:]
  for i in range(ncycles):
    times = np.append(times, times[-1] + np.linspace(0, per, nsteps_cycle+1)[1:])
  
  times = np.insert(times, 0, 0)
  dts = np.diff(times)

  model.add_force_bc(2, ffn)
  model.add_displacement_bc(1, fixed)
  
  for dt in dts:
    model.solve(dt)
  
  return model.node[2]['displacements'][-1]
  
if __name__ == "__main__":
  E = 200000.0
  nu = 0.27
  sy = 300.0
  k = E / 100
  a = 1.0e-5

  mat = make_mat(E, nu, sy, k, a)

  A1 = 10.0
  A2 = 17.0

  l = 125.0

  P = 2500.0
  dT = 200.0
  ncycles = 20
  
  d_old = old_way(A1, A2, mat, l, P, dT, ncycles)
  d_new = new_way(A1, A2, mat, l, P, dT, ncycles)

  print("Old way: %f" % d_old)
  print("New way: %f" % d_new)
