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
def simplify_Ti_model(x_sample, taus_1, taus_2, taus_3,
                      taut_1, taut_2, X_s, 
                      k1_1, k1_2, k1_3, X, 
                      k2_1, k2_2, k2_3 
                      ):
#================================================#
  
  res = simplify_model(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s,
            k1_1, k1_2, k1_3, X, 
            k2_1, k2_2, k2_3,
            T = 298.0, emax = emax, N = 10, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True, PTR = True)
            
  return res

#================================================#
def make_Ti_model(x_sample, taus_1, taus_2, taus_3,
                  taut_1, taut_2, X_s, 
                  k1_1, k1_2, k1_3, 
                  X, g_1, g_2, g_3, tau_D1, tau_D2, 
                  tau_D3):
#================================================#
  
  res = make_model(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s, 
            k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = emax, N = 5, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True, Taylor = True,
            PTR = True)

  return res



#================================================#
def interpolate_obs():
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1) 
  return df

  
if __name__ == "__main__":
  
  simple = False

  if simple:
    params = [270.0, 190.0, 310.0,
            280.0, 350.0, 0.9,
            2.00, 2.00, 5.00,
            0.35, 70.0, 70.0, 70.0]
    df = interpolate_obs()
    res = simplify_Ti_model(x_sample, *params)
    plt.plot(res['strain'], res['stress'], 'g--', label = 'fit')
    plt.plot(df['Nominal_strain'], df['True_stress'], 'b-', label = 'data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.savefig("Ti-fitting-simple.png")
    plt.show()
    plt.close()
    
  else:
    params = [270.0, 190.0, 310.0,
            280.0, 350.0, 0.9,
            2.00, 2.00, 5.00,
            0.35, 0.0045, 0.0045, 0.005,
            100.0, 100.0, 100.0]
    df = interpolate_obs()
    res = make_Ti_model(x_sample, *params)
    plt.plot(res['strain'], res['stress'], 'g--', label = 'fit')
    plt.plot(df['Nominal_strain'], df['True_stress'], 'b-', label = 'data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.savefig("Ti-fitting-complex.png")
    plt.show()
    plt.close()
  
