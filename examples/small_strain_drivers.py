#!/usr/bin/env python3

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow, visco_flow, general_flow, interpolate

import matplotlib.pyplot as plt
import numpy as np

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

  # Uniaxial tension 
  erate = 1.0e-4
  res = drivers.uniaxial_test(model, erate)
  plt.figure()
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.title("Uniaxial tension")
  plt.show()
  
  # Strain controlled-cyclic
  emax = 0.01
  R = -0.75
  erate = 1.0e-4
  ncycles = 5
  hold_time = [10*3600.0,0]
  res = drivers.strain_cyclic(model, emax, R, erate, ncycles, 
      hold_time = hold_time)
  plt.figure()
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.title("Strain-controlled cyclic")
  plt.show()

  # Strain controlled-cyclic with follow up
  emax = 0.01
  R = -0.75
  erate = 1.0e-4
  ncycles = 5
  hold_time = [10*3600.0,0]
  res = drivers.strain_cyclic_followup(model, emax, R, erate, ncycles, 
      hold_time = hold_time, q = 4)
  plt.figure()
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.title("Strain-controlled cyclic with follow up")
  plt.show()

  # Stress controlled-cyclic
  smax = 200.0
  R = -0.75
  erate = 1.0e-4
  ncycles = 5
  hold_time = [10*3600.0,0]
  res = drivers.stress_cyclic(model, smax, R, erate, ncycles, 
      hold_time = hold_time)
  plt.figure()
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.title("Stress-controlled cyclic")
  plt.show()

  # Stress relaxation
  emax = 0.05
  erate = 1.0e-4
  hold = 1000.0*3600.0
  res = drivers.stress_relaxation(model, emax, erate, hold)
  plt.figure()
  plt.plot(res['rtime']/3600.0, res['rstress'], 'k-')
  plt.xlabel("Time (hrs)")
  plt.ylabel("Stress (MPa)")
  plt.title("Stress relaxation")
  plt.tight_layout()
  plt.show()

  # Creep
  smax = 180.0
  srate = 1.0
  hold = 10000.0 * 3600.0
  res = drivers.creep(model, smax, srate, hold, logspace = True)
  plt.figure()
  plt.plot(res['rtime']/3600.0, res['rstrain'], 'k-')
  plt.xlabel("Time (hrs)")
  plt.ylabel("Creep strain (mm/mm)")
  plt.title("Creep")
  plt.tight_layout()
  plt.show()

  # Rate jump
  erates = [1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e-2, 1.0e-4, 1.0e-6]
  res = drivers.rate_jump_test(model, erates)
  plt.figure()
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.title("Strain rate jump test")
  plt.tight_layout()
  plt.show()

  # Isochronous curve
  time = 50000.0 * 3600.0
  res = drivers.isochronous_curve(model, time)
  plt.figure()
  plt.plot(res['strain'], res['stress'], 'k-')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.title("Isochronous stress-strain curve")
  plt.tight_layout()
  plt.show()
  
  # Need a temperature sensitive model for this one!
  vmodel = visco_flow.YaguchiGr91FlowRule()

  E_poly = np.array([-3.077e-1, 3.003e2, 1.269e5])
  nu = 0.3

  mu_poly = list(E_poly/(2*(1+nu)))
  K_poly = list(E_poly/(3*(1-2*nu)))

  shear = interpolate.PolynomialInterpolate(mu_poly)
  bulk = interpolate.PolynomialInterpolate(K_poly)
  elastic = elasticity.IsotropicLinearElasticModel(shear, "shear", 
      bulk, "bulk")

  flow = general_flow.TVPFlowRule(elastic, vmodel)

  model = models.GeneralIntegrator(elastic, flow)
  
  # Thermomechanical
  time = list(np.linspace(0, 1000))
  Ts = list(np.linspace(500+273.15,600+273.15))
  strains = list(np.linspace(0.0,0.02))
  
  time = time + time[::-1][1:]
  Ts = Ts + Ts[::-1][1:]
  strains = strains + strains[::-1][1:]

  res = drivers.thermomechanical_strain_raw(model, time, Ts,
      strains)

  fig, ax1 = plt.subplots()
  ax1.plot(res['strain'], res['stress'], 'k-', label = "stress")
  ax1.set_xlabel("Strain (mm/mm)")
  ax1.set_ylabel("Stress (MPa)")
  
  ax2 = ax1.twinx()
  ax2.plot(res['strain'], res['temperature'], 'r-', label = "temperature")
  ax2.set_ylabel("Temperature (K)")

  ax1.legend(loc = 'best')
  ax2.legend(loc = 'best')

  plt.title("Thermomechanical test")
  plt.tight_layout()
  plt.show()
