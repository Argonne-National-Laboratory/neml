#!/usr/bin/env python

import numpy as np
import scipy.interpolate as inter
import scipy.optimize as opt
import os.path

import matplotlib.pyplot as plt

import sys
sys.path.append('../..')

from neml import solvers, neml, interpolate, elasticity, drivers, surfaces, hardening, ri_flow, visco_flow, general_flow

from generate import *

if __name__ == "__main__":
  temperature_C = np.array([500, 550])

  stress = 200.0 # MPa
  stress_rate = 1.0e-4 # (MPa/s)
  hold = 50000.0 * 3600.0 # 500000 hours in seconds

  # Going to plot 2 different temperatures as blue and red lines
  # Make the Yaguchi model solid and the Koo model dashed lines
  colors = ['b', 'r']
  
  # Yaguchi model interpolates over temperature, only need 1 instantiation
  yaguchi = make_yaguchi_model()

  for tC, color in zip(temperature_C, colors):
    tK = tC + 273.15
    # For the Koo model call a function to make an instantiation of a 
    # Chaboche model with the appropriate temperature-dependent constants
    koo = make_koo_model(str(tC))

    # Run the yaguchi model
    resy = drivers.creep(yaguchi, stress, stress_rate, hold, T = tK)

    # Run the Koo model
    resk = drivers.creep(koo, stress, stress_rate, hold) # T doesn't matter
    
    # Plot the results for this temperature
    plt.loglog(resy['rtime'] / 3600.0, resy['rrate'], color+'-')
    plt.loglog(resk['rtime'] / 3600.0, resk['rrate'], color+'--')

  plt.xlabel("Time (hrs)")
  plt.ylabel("Creep rate (1/s)")
  plt.show()
