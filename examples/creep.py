#!/usr/bin/env python

import sys
sys.path.append('..')

import numpy as np
import scipy.integrate as quad

from neml import solvers, creep

import matplotlib.pyplot as plt

if __name__ == "__main__":
  # T is in hours, strain in percent, stress in MPa
  A = 1.85e-9
  n = 2.5
  m = 0.3

  smodel = creep.NortonBaileyCreep(A, n, m)
  model = creep.J2CreepModel(smodel)

  stress = 150.0
  temp = 300.0
  times = np.linspace(0, 500, 100)

  # Function to integrate out the scalar model
  def scalar_rate(e, t):
    return smodel.g(stress, e[0], t, temp)
  
  estart = 1.0e-14
  strains = quad.odeint(scalar_rate, estart, times)
  plt.plot(times, strains, 'k-')

  vstrains = []
  vstresses = []

  stress = np.array([stress, 0, 0, 0, 0, 0])
  vstresses.append(stress)
  vstrains.append(np.array([estart, -0.5*estart, -0.5*estart, 0, 0, 0]))
  for i in range(1, len(times)):
    strain_next, A_next = model.update(stress, vstrains[-1], temp, temp, 
        times[i], times[i-1])
    vstresses.append(stress[-1])
    vstrains.append(strain_next)

  plt.show()
