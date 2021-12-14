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


# sets up x_scale for both experiment and simulation
emax = 0.2
Nsample = 200
x_sample = np.linspace(0.0, emax*0.99, Nsample)



#================================================#
def simplify_Ti_model(x_sample, X_s, k1_1, k1_2, k1_3, X, 
                      k2_1, k2_2, k2_3, 
                      check = True):
#================================================#
  
  res = simplify_model(X_s, k1_1, k1_2, k1_3, X, 
            k2_1, k2_2, k2_3,
            T = 298.0, emax = emax, N = 5, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True, PTR = True)
  yobs = interpolate(res['strain'], res['stress'], x_sample)
  if check:
    return yobs
  else:
    return res

#================================================#
def make_Ti_model(x_sample, X_s, k1_1, k1_2, k1_3, 
                  X, g_1, g_2, g_3, tau_D1, tau_D2, 
                  tau_D3, check = True):
#================================================#
  
  res = make_model(X_s, k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = emax, N = 5, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True,
            PTR = True)
  yobs = interpolate(res['strain'], res['stress'], x_sample)
  if check:
    return yobs
  else:
    return res



#================================================#
def interpolate_obs():
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1) 
  return df

  
if __name__ == "__main__":
  
  check = False
  fa, fac = 1.0, 1.0

  params = [0.9, 1.00000009*fa, 0.99999987*fa, 
            0.99999987*fa, 0.5, 1.00000002*fac,
            1.00000013*fac, 1.00000001*fac]
  """
  params = [0.99999998, 3.0, 3.0,
            3.0, 0.99999988, 0.016, 0.016, 0.016, 
            100.0, 100.0, 100.0]
  """
  if check:
    df = interpolate_obs()
    y_sim = make_Ti_model(x_sample, *params, check = check)
    plt.plot(x_sample, y_sim, 'g--', label = 'fit')
    plt.plot(df['Nominal_strain'], df['True_stress'], 'b-', label = 'data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.savefig("Ti-fitting-manual.png")
    plt.show()
    plt.close()
    
  else:
    df = interpolate_obs()
    res = simplify_Ti_model(x_sample, *params, check = check)
    plt.plot(res['strain'], res['stress'], 'g--', label = 'fit')
    plt.plot(df['Nominal_strain'], df['True_stress'], 'b-', label = 'data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()
    plt.close()
  
