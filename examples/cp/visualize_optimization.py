#!/usr/bin/env python3

import sys
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


# sets up x_scale for both experiment and simulation
emax = 0.2
Nsample = 200
x_sample = np.linspace(0.0, emax*0.99, Nsample)

# sets up the parameters range 
min_theta = (X_s_min, k1_1_min, k1_2_min, k1_3_min, 
        X_min, g_1_min, g_2_min, g_3_min,
        tau_D1_min, tau_D2_min, tau_D3_min) = (3.150832905, 41.90999565,
                                               9.04691245, 122.1456805,
                                               0.4906251715, 0.0145345937,
                                               0.0114694439, 0.03620099355,
                                               53.5702925, 50.0646285,
                                               45.9480411)
                                               
max_theta = (X_s_max, k1_1_max, k1_2_max, k1_3_max, 
        X_max, g_1_max, g_2_max, g_3_max,
        tau_D1_max, tau_D2_max, tau_D3_max) = (4.7262493575, 125.72998695,
                                               27.140737350000002, 366.43704149999996,
                                               1.4718755145, 0.0436037811,
                                               0.034408331699999996, 0.10860298064999999,
                                               160.7108775, 150.1938855,
                                               137.84412329999998)
"""

# sets up initial model parameters
(X_s_i, k1_1_i, k1_2_i, k1_3_i, 
        X_i, g_1_i, g_2_i, g_3_i, 
        tau_D1_i, tau_D2_i, tau_D3_i) = (0.9, 75.0, 12.5, 250.0,
                                         0.3, 0.02, 0.02, 0.055,
                                         100.0, 100.0, 100.0)
                                         
sf = 0.99

# sets up x_scale for both experiment and simulation
emax = 0.05
Nsample = 200
x_sample = np.linspace(0.0, emax*0.99, Nsample)

# sets up the parameters range 
min_theta = (X_s_min, k1_1_min, k1_2_min, k1_3_min, 
        X_min, g_1_min, g_2_min, g_3_min,
        tau_D1_min, tau_D2_min, tau_D3_min) = (X_s_i*(1-sf), k1_1_i*(1-sf), 
                                               k1_2_i*(1-sf), k1_3_i*(1-sf), 
                                               X_i*(1-sf), g_1_i*(1-sf),
                                               g_2_i*(1-sf), g_3_i*(1-sf), 
                                               tau_D1_i*(1-sf), tau_D2_i*(1-sf),
                                               tau_D3_i*(1-sf))
                                               
max_theta = (X_s_max, k1_1_max, k1_2_max, k1_3_max, 
        X_max, g_1_max, g_2_max, g_3_max,
        tau_D1_max, tau_D2_max, tau_D3_max) = (X_s_i*(1+sf), k1_1_i*(1+sf), 
                                               k1_2_i*(1+sf), k1_3_i*(1+sf), 
                                               X_i*(1+sf), g_1_i*(1+sf),
                                               g_2_i*(1+sf), g_3_i*(1+sf), 
                                               tau_D1_i*(1+sf), tau_D2_i*(1+sf),
                                               tau_D3_i*(1+sf))
"""
#================================================#
def convert_to_real(p):
#================================================#
    bounds = np.array([
                    [X_s_min, X_s_max],    # X_s
                    [k1_1_min, k1_1_max],    # k1_1
                    [k1_2_min, k1_2_max],  # k1_2
                    [k1_3_min, k1_3_max],  # k1_3
                    [X_min, X_max],  # X
                    [g_1_min, g_1_max],  # g_1
                    [g_2_min, g_2_max],  # g_2
                    [g_3_min, g_3_max],  # g_3
                    [tau_D1_min, tau_D1_max],  # tau_D1
                    [tau_D2_min, tau_D2_max],  # tau_D2
                    [tau_D3_min, tau_D3_max],  # tau_D3
                 ])

    return bounds[:,0] + (p * (bounds[:,1] - bounds[:,0]))

#================================================#
def make_Ti_model(params):
#================================================#
  
  theta_in = (params[0], 
         params[1], params[2], params[3], 
         params[4], params[5], params[6],
         params[7], params[8], params[9],
         params[10])
  
  theta = convert_to_real(theta_in)
  
  X_s, k1_1, k1_2, k1_3, X, g_1, g_2, g_3, tau_D1, tau_D2, tau_D3  = (theta[0],theta[1],
        theta[2],theta[3],theta[4],theta[5],theta[6],theta[7],theta[8],theta[9],theta[10])
  
  res = make_model(X_s, k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = emax, N = 10, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = False, Taylor = True, 
            PTR = False)

  return res


#================================================#
def interpolate_obs(x_sample):
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1)
  stress = interpolate(df['Nominal_strain'], df['True_stress'], x_sample)  
  return stress




if __name__ == "__main__":
  
  T = 298.0

  params = [0.28776328, 0.48671358, 0.74433297, 0.2457169,
            0.07744721, 0.92449378, 0.17374999, 0.28382545,
            0.02952822, 0.10566986, 0.65129651]
  theta_in = (params[0], 
         params[1], params[2], params[3], 
         params[4], params[5], params[6],
         params[7], params[8], params[9],
         params[10])
  
  res = make_Ti_model(params)
  stress = interpolate_obs(x_sample)
  res_stress = interpolate(res['strain'], res['stress'], x_sample)
  
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1)
  
  plt.plot(res['strain'], res['stress'], label = "Sim - %3.0f C" % (T-273.15))
  # plt.plot(x_sample, res_stress, label = "Sim - %3.0f C interpolate" % (T-273.15))
  plt.plot(df['Nominal_strain'], df['True_stress'], label = "Exp - %3.0f C" % (T-273.15))
  # plt.plot(x_sample, stress, label = "Exp - %3.0f C interpolate" % (T-273.15))

  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  # plt.savefig("tension-Ti.png")
  plt.show()
  plt.close()
