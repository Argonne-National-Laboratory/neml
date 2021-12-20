#!/usr/bin/env python3

from Ti_model_make import *
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


# sets up x_scale for both experiment and simulation
emax = 0.4
Nsample = 1000
x_sample = np.linspace(0.0, emax*0.99, Nsample)

#  model grains and threads
erate = 1.0e-3
Ngrains = 200
nthreads = 30
T = 373.0


#================================================#
def fit_Ti_temp(x_sample, taus_1, taus_2, taus_3,
                taut_1, taut_2, k2_1, k2_2, k2_3):
#================================================#
  
  res = Ti_maker_sim(taus_1, taus_2, taus_3,
            taut_1, taut_2, 0.9,
            1.00, 0.15, 4.00,
            k2_1, k2_2, k2_3,
            T = T, emax = emax, N = Ngrains, 
            strain_rate = erate, nthreads = nthreads, 
            verbose = True, Taylor = True, PTR = True)
            
  yobs = interp(res['strain'], res['stress'], x_sample)
  return yobs


#================================================#
def fit_Ti_model(x_sample, taus_1, taus_2, taus_3,
                         taut_1, taut_2, X_s, 
                         k1_1, k1_2, k1_3, 
                         X, k2_1, k2_2, k2_3):
#================================================#
  
  res = Ti_maker_sim(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s,
            k1_1, k1_2, k1_3,
            k2_1, k2_2, k2_3,
            T = T, emax = emax, N = Ngrains, 
            strain_rate = erate, nthreads = nthreads, 
            verbose = True, Taylor = True, PTR = True)
            
  yobs = interp(res['strain'], res['stress'], x_sample)
  return yobs


#================================================#
def load_data(T):
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Nasser-1999-Acta/"
  df = load_file(path_1, str(int(T))) 
  return df

  
if __name__ == "__main__":
  
  # sets up parameters range
  min_params = [95.0, 75.5, 95.0, 120.0, 190.0, 16.0, 16.0, 16.0]
  max_params = [100.0, 80.5, 100.0, 130.0, 200.0, 25.0, 25.0, 25.0]

  df = load_data(T)
  xdata = x_sample
  ydata = interp(df['Nominal_strain'], df['True_stress'], x_sample)
  popt, pcov = curve_fit(fit_Ti_temp, xdata, ydata)
    bounds=(min_params, max_params))

  print("popt: ", popt)
  print("")
  print("pcov: ", pcov)

  plt.plot(xdata, fit_Ti_temp(xdata, *popt), 'g--', label = 'fit')
     # label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
  plt.plot(xdata, ydata, 'b-', label='data')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.legend()
  plt.grid(True)
  # plt.savefig("tension-Ti-simplify.png")
  plt.show()
  plt.close()

  