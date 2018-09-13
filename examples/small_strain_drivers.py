#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow, visco_flow, general_flow

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

def example_econt_erate(model, sdir, rate, T, etotal, nsteps):
  """
    Parameters:
      model     material model
      sdir      stress direction
      rate      strain rate
      T         constant temperature
      etotal    total strain, in direction
      nsteps    number of steps
  """
  e_inc = etotal / nsteps

  driver = drivers.Driver_sd(model, verbose = True)

  print("Rate controlled...")
  for i in range(nsteps):
    print(i+1)
    driver.erate_einc_step(sdir, rate, e_inc, T)
  
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

def example_scont_srate(model, sdir, rate, T, stotal, nsteps):
  """
    Parameters:
      model     material model
      sdir      stress direction
      rate      stress rate
      T         constant temperature
      etotal    total stress, in direction
      nsteps    number of steps
  """
  s_inc = stotal / nsteps

  driver = drivers.Driver_sd(model, verbose = True)

  print("Rate controlled...")
  for i in range(nsteps):
    print(i+1)
    driver.srate_sinc_step(sdir, rate, s_inc, T)
  
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
  n = 20.0
  eta = 108.0
  sY = 89.0

  Q = 165.0
  b = 12.0
  
  C1 = 80.0e3
  C2 = 14.02e3
  C3 = 3.333e3

  y1 = 0.9e3
  y2 = 1.5e3
  y3 = 1.0

  surface = surfaces.IsoKinJ2()
  iso = hardening.VoceIsotropicHardeningRule(sY, Q, b)
  cs = [C1, C2, C3]
  gs = np.array([y1, y2, y3])
  gmodels = [hardening.ConstantGamma(g) for g in gs]
  hmodel = hardening.Chaboche(iso, cs, gmodels, [0.0] * len(cs),
      [1.0] * len(cs))

  fluidity = visco_flow.ConstantFluidity(eta)

  vmodel = visco_flow.ChabocheFlowRule(surface, hmodel, fluidity, n)

  E = 92000.0
  nu = 0.3

  mu = E/(2*(1+nu))
  K = E/(3*(1-2*nu))

  elastic = elasticity.IsotropicLinearElasticModel(mu, "shear", K, 
      "bulk")

  flow = general_flow.TVPFlowRule(elastic, vmodel)

  model = models.GeneralIntegrator(elastic, flow, verbose = False) 

  example_strain(model, np.array([0.04,-0.02,-0.02,0,0,0]), 300.0, 100.0, 100)
  example_stress(model, np.array([300.0,0,0,0,0,0]), 300.0, 100, 100)
  example_econt_erate(model, np.array([1,0,0,0,0,0]), 1.0e0, 300.0, 0.04, 100)
  example_scont_srate(model, np.array([1,0,0,0,0,0]), 1.0, 300.0, 300.0, 100)
