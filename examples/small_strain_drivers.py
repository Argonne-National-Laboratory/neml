#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import neml, elasticity, drivers, surfaces, hardening, ri_flow

import matplotlib.pyplot as plt
import numpy as np

def example_strain(model, strain, T, t, nsteps):
  """
    Parameters:
      model     material model
      strain    target strain
      T         constant temperature
      t         time to take
      nsteps    number of steps
  """
  de = strain / nsteps
  dt = t / nsteps
  
  e = np.zeros((6,))
  t = 0

  driver = drivers.Driver_sd(model, verbose = True)
  
  print("Driving strain...")
  for i in range(nsteps):
    print(i+1)
    e += de
    t += dt

    driver.strain_step(np.copy(e), t, T)

  plt.plot(driver.strain[:,0], driver.stress[:,0], 'k-')
  plt.plot(driver.strain[:,0], driver.stress[:,1], 'r-')
  plt.plot(driver.strain[:,0], driver.stress[:,2], 'b-')

  plt.plot(driver.strain[:,0], driver.stress[:,3], 'k--')
  plt.plot(driver.strain[:,0], driver.stress[:,4], 'r--')
  plt.plot(driver.strain[:,0], driver.stress[:,5], 'b--')
  plt.show()

  plt.plot(driver.strain[:,0], driver.history[:,0], 'k-')
  plt.plot(driver.strain[:,0], driver.history[:,1], 'r-')
  plt.plot(driver.strain[:,0], driver.history[:,2], 'b-')

  plt.plot(driver.strain[:,0], driver.history[:,3], 'k--')
  plt.plot(driver.strain[:,0], driver.history[:,4], 'r--')
  plt.plot(driver.strain[:,0], driver.history[:,5], 'b--')
  plt.show()

  plt.plot(driver.strain[:,0], driver.history[:,6], 'k-')
  plt.show()

def example_stress(model, stress, T, t, nsteps):
  """
    Parameters:
      model     material model
      stress    target stress
      T         constant temperature
      t         time to take
      nsteps    number of steps
  """
  ds = stress / nsteps
  dt = t / nsteps
  
  s = np.zeros((6,))
  t = 0

  driver = drivers.Driver_sd(model, verbose = True)
  
  print("Driving stress...")
  for i in range(nsteps):
    print(i+1)
    s += ds
    t += dt

    driver.stress_step(s, t, T)
 
  plt.plot(driver.strain[:,0], driver.stress[:,0], 'k-')
  plt.plot(driver.strain[:,0], driver.stress[:,1], 'r-')
  plt.plot(driver.strain[:,0], driver.stress[:,2], 'b-')

  plt.plot(driver.strain[:,0], driver.stress[:,3], 'k--')
  plt.plot(driver.strain[:,0], driver.stress[:,4], 'r--')
  plt.plot(driver.strain[:,0], driver.stress[:,5], 'b--')
  plt.show()

  plt.plot(driver.strain[:,0], driver.history[:,0], 'k-')
  plt.plot(driver.strain[:,0], driver.history[:,1], 'r-')
  plt.plot(driver.strain[:,0], driver.history[:,2], 'b-')

  plt.plot(driver.strain[:,0], driver.history[:,3], 'k--')
  plt.plot(driver.strain[:,0], driver.history[:,4], 'r--')
  plt.plot(driver.strain[:,0], driver.history[:,5], 'b--')
  plt.show()

  plt.plot(driver.strain[:,0], driver.history[:,6], 'k-')
  plt.show()

def example_rate(model, sdir, rate, T, dt, nsteps):
  """
    Parameters:
      model     material model
      sdir      stress direction
      rate      strain rate
      T         constant temperature
      dt        time increment
      nsteps    number of steps
  """
  t = 0

  driver = drivers.Driver_sd(model, verbose = True)

  print("Rate controlled...")
  for i in range(nsteps):
    print(i+1)
    t += dt
    driver.rate_step(sdir, rate, t, T)
  
  plt.plot(driver.strain[:,0], driver.stress[:,0], 'k-')
  plt.plot(driver.strain[:,0], driver.stress[:,1], 'r-')
  plt.plot(driver.strain[:,0], driver.stress[:,2], 'b-')

  plt.plot(driver.strain[:,0], driver.stress[:,3], 'k--')
  plt.plot(driver.strain[:,0], driver.stress[:,4], 'r--')
  plt.plot(driver.strain[:,0], driver.stress[:,5], 'b--')
  plt.show()

  plt.plot(driver.strain[:,0], driver.history[:,0], 'k-')
  plt.plot(driver.strain[:,0], driver.history[:,1], 'r-')
  plt.plot(driver.strain[:,0], driver.history[:,2], 'b-')

  plt.plot(driver.strain[:,0], driver.history[:,3], 'k--')
  plt.plot(driver.strain[:,0], driver.history[:,4], 'r--')
  plt.plot(driver.strain[:,0], driver.history[:,5], 'b--')
  plt.show()

  plt.plot(driver.strain[:,0], driver.history[:,6], 'k-')
  plt.show()

if __name__ == "__main__":
  E = 200000.0
  nu = 0.3

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  s0 = 150.0
  R = 100.0
  d = 1000.0
  
  shear = elasticity.ConstantShearModulus(mu)
  bulk = elasticity.ConstantBulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)
  #model = neml.SmallStrainElasticity(elastic)
  surface = surfaces.IsoJ2()
  #hrule = hardening.LinearIsotropicHardeningRule(s0, Kp)
  hrule = hardening.VoceIsotropicHardeningRule(s0, R, d)
  flow = ri_flow.RateIndependentAssociativeFlow(surface, hrule)
  model = neml.SmallStrainRateIndependentPlasticity(elastic, flow)

  example_strain(model, np.array([0.01,0,0,0,0,0]), 300.0, 10, 100)
  example_stress(model, np.array([220.0,0,0,0,0,0]), 300.0, 10, 20)
  example_rate(model, np.array([1,0,0,0,0,0]), 1.0e-2, 300.0, 2.0e-3, 100)

