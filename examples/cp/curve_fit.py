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
emax = 0.15
Nsample = 1000
x_sample = np.linspace(0.0, emax*0.99, Nsample)


#================================================#
def make_Ti_simple_model(x_sample, taus_1, taus_2, taus_3,
                         taut_1, taut_2, X_s, 
                         k1_1, k1_2, k1_3, 
                         X, k2_1, k2_2, k2_3):
#================================================#
  
  res = simplify_model(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s,
            k1_1, k1_2, k1_3, X, 
            k2_1, k2_2, k2_3,
            T = 298.0, emax = emax, N = 20, 
            strain_rate = 1.0e-4, nthreads = 30, 
            verbose = True, Taylor = True, PTR = True)
            
  yobs = interpolate(res['strain'], res['stress'], x_sample)
  return yobs

#================================================#
def make_Ti_model(x_sample, taus_1, taus_2, taus_3,
                  taut_1, taut_2, X_s, 
                  k1_1, k1_2, k1_3, X, 
                  g_1, g_2, g_3,
                  tau_D1, tau_D2, tau_D3):
#================================================#
  
  res = make_model(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s, 
            k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = emax, N = 20, 
            strain_rate = 1.0e-4, nthreads = 30, 
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
  
  simplify = False
  
  if simplify:
    # sets up parameters range
    min_params = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    max_params = [1.0, 100.0, 100.0, 100.0, 1.0, 100.0, 100.0, 100.0]

    df = interpolate_obs()
    xdata = x_sample
    ydata = interpolate(df['Nominal_strain'], df['True_stress'], x_sample)
    popt, pcov = curve_fit(make_Ti_simple_model, xdata, ydata)
      # bounds=(min_params, max_params))

    print("popt: ", popt)
    print("")
    print("pcov: ", pcov)

    plt.plot(xdata, make_Ti_simple_model(xdata, *popt), 'g--', label = 'fit')
         # label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.savefig("tension-Ti-simplify.png")
    plt.show()
    plt.close()

  else:
  
    # sets up parameters range
    min_params = [0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 100.0, 100.0, 100.0]
    max_params = [1.0, 100.0, 100.0, 100.0, 1.0, 1.0, 1.0, 1.0, 10000.0, 10000.0, 10000.0]

    df = interpolate_obs()
    xdata = x_sample
    ydata = interpolate(df['Nominal_strain'], df['True_stress'], x_sample)
    popt, pcov = curve_fit(make_Ti_model, xdata, ydata)
      # bounds=(min_params, max_params))

    print("popt: ", popt)
    print("")
    print("pcov: ", pcov)

    plt.plot(xdata, make_Ti_model(xdata, *popt), 'g--', label = 'fit')
         # label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    plt.plot(xdata, ydata, 'b-', label='data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.savefig("tension-Ti-allparam.png")
    plt.show()
    plt.close()
  