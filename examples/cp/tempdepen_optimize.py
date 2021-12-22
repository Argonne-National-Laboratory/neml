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
def fit_Ti_temp(taus_1, taus_2, taus_3,
                taut_1, taut_2, k2_1, k2_2, k2_3):
#================================================#
  
  res = Ti_maker_sim(taus_1, taus_2, taus_3,
            taut_1, taut_2, 0.9,
            1.00, 0.15, 4.00,
            k2_1, k2_2, k2_3,
            T = T, emax = emax, N = Ngrains, 
            strain_rate = erate, nthreads = nthreads, 
            verbose = True, Taylor = True, PTR = True)
            
  return res

#================================================#
def load_data(T):
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Nasser-1999-Acta/"
  df = load_file(path_1, str(int(T))) 
  return df


#================================================#
def R(params):
#================================================#
  res = fit_Ti_temp(*params)
  df = load_data(T)
  yobs = interp(df['Nominal_strain'], df['True_stress'], x_sample)
  ysim = interp(res['strain'], res['stress'], x_sample)
  R = la.norm(ysim - yobs)
  print("Current residual: %e" % R)
  return R
  
if __name__ == "__main__":
  
  # sets up the range of model parameters
  taus_1_range = [95.0, 100.0]
  taus_2_range = [75.0, 80.5]
  taus_3_range = [95.0, 100.0]
  taut_1_range = [120.0, 130.0]
  taut_2_range = [190.0, 200.0]
  k2_1_range = [16.0, 25.0]
  k2_2_range = [16.0, 25.0]
  k2_3_range = [16.0, 25.0]

  p0 = [ra.uniform(*taus_1_range), ra.uniform(*taus_2_range),
        ra.uniform(*taus_3_range), ra.uniform(*taut_1_range),
        ra.uniform(*taut_2_range), ra.uniform(*k2_1_range),
        ra.uniform(*k2_2_range), ra.uniform(*k2_3_range)]

  # res = minimize_parallel(R, p0, bounds = [taus_1_range,
  res = opt.minimize(R, p0, method = 'L-BFGS-B', bounds = [taus_1_range,
        taus_2_range, taus_3_range, taut_1_range, taut_2_range,
        k2_1_range, k2_2_range, k2_3_range])

  print(res.success)
  print(res)

  ref_n = (taus_1_n, taus_2_n, taus_3_n, taut_1_n, taut_2_n, k2_1_n, k2_2_n, k2_3_n) = (res.x[0],
                res.x[1], res.x[2], res.x[3], res.x[4],
                res.x[5], res.x[6], res.x[7])
                
  res_final = fit_Ti_temp(*ref_n)
  df = load_data(T)
  # visualize the fitting
  plt.plot(res_final['strain'], res_final['stress'], label = "Sim - %3.0f C" % (T-273.15))
  plt.plot(df['Nominal_strain'], df['True_stress'], label = "Exp - %3.0f C" % (T-273.15))
  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  # plt.savefig("tension-Ti.png")
  plt.show()
  plt.close()