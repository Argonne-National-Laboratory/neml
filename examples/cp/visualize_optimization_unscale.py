#!/usr/bin/env python3

from Ti_uniaxial_fitting import *
from neml import drivers

import numpy as np
import numpy.linalg as la
import numpy.random as ra
import scipy.interpolate as inter
import scipy.optimize as opt
import matplotlib.pyplot as plt
import pandas as pd

import time
import concurrent.futures
from multiprocessing import Pool
from optimparallel import minimize_parallel

from scipy.optimize import curve_fit


sf = 0.8

# sets up x_scale for both experiment and simulation
emax = 0.2
Nsample = 200
x_sample = np.linspace(0.0, emax*0.99, Nsample)


#================================================#
def make_Ti_model(x_sample, X_s, k1_1, k1_2, k1_3, 
                  X, g_1, g_2, g_3, tau_D1, tau_D2, 
                  tau_D3):
#================================================#
  
  res = make_model(X_s, k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = emax, N = 50, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True,
            PTR = True)
  yobs = interpolate(res['strain'], res['stress'], x_sample)
  return yobs

#================================================#
def interpolate_obs():
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1) 
  return df

  
if __name__ == "__main__":
  
  df = interpolate_obs()
  xdata = x_sample
  ydata = interpolate(df['Nominal_strain'], df['True_stress'], x_sample)
  params = [1.00000176, 0.99999989, 0.99999972, 1.00000035, 0.99999971,
            1.00000174, 0.99999919, 0.99999946, 0.9999988, 1.00000011,
            0.99999834]
  y_sim = make_Ti_model(x_sample, *params)
  
  
  plt.plot(x_sample, y_sim, 'g--', label = 'fit')
  plt.plot(xdata, ydata, 'b-', label='data')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.legend()
  plt.show()
  