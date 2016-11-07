#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import neml, volumetric, deviatoric, shear, drivers, surface, hardening

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

    driver.strain_step(e, t, T)

  plt.plot(driver.strain[:,0], driver.stress[:,0], 'k-')
  plt.plot(driver.strain[:,0], driver.stress[:,1], 'r-')
  plt.plot(driver.strain[:,0], driver.stress[:,2], 'b-')

  plt.plot(driver.strain[:,0], driver.stress[:,3], 'k--')
  plt.plot(driver.strain[:,0], driver.stress[:,4], 'r--')
  plt.plot(driver.strain[:,0], driver.stress[:,5], 'b--')
  plt.show()
 
  plt.plot(driver.strain[:,0], driver.history[:,0], 'k-')
  plt.show()

  plt.plot(driver.strain[:,0], driver.history[:,1], 'k-')
  plt.plot(driver.strain[:,0], driver.history[:,2], 'r-')
  plt.plot(driver.strain[:,0], driver.history[:,3], 'b-')

  plt.plot(driver.strain[:,0], driver.history[:,4], 'k--')
  plt.plot(driver.strain[:,0], driver.history[:,5], 'r--')
  plt.plot(driver.strain[:,0], driver.history[:,6], 'b--')
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

if __name__ == "__main__":
  E = 200000.0
  nu = 0.3

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))
  
  K0 = 200.0
  Kb = E / 100.0
  #Ksat = 250.0
  #delta = 100.0

  vol_model = volumetric.VModelConstantK(K)
  shear_model = shear.ConstantShearModulus(mu)
  ys = surface.IsoJ2()
  hr = hardening.IsoJ2LinearAHardening(K0, Kb)
  dev_model = deviatoric.RIAFModel(shear_model, ys, hr, verbose = False)
  #dev_model = deviatoric.LEModel(shear_model)
  model = neml.SplitModel_sd(vol_model, dev_model)

  example_strain(model, np.array([0.02,0,0,0,0,0]), 300.0, 10, 100)
  example_stress(model, np.array([220.0,0,0,0,0,0]), 300.0, 10, 20)
  #example_rate(model, np.array([1,0,0,0,0,0]), 1.0e-2, 300.0, 1.0e-1, 50)

