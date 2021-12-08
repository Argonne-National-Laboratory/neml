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

"""
# sets up x_scale for both experiment and simulation
emax = 0.05
Nsample = 200
x_sample = np.linspace(0.0, emax*0.99, Nsample)

# sets up the parameters range 
min_theta = (X_s_min, k1_1_min, k1_2_min, k1_3_min, 
        X_min, g_1_min, g_2_min, g_3_min,
        tau_D1_min, tau_D2_min, tau_D3_min) = (0.527777705, 74.625,
                                               5.30971905, 61.379739,
                                               0.255167712, 0.01367479605,
                                               0.01396919455, 0.01855309345,
                                               49.6123269, 90.436568,
                                               42.58838775)
                                               
max_theta = (X_s_max, k1_1_max, k1_2_max, k1_3_max, 
        X_max, g_1_max, g_2_max, g_3_max,
        tau_D1_max, tau_D2_max, tau_D3_max) = (1.5833331149999998, 223.875,
                                               15.92915715, 184.139217,
                                               0.765503136, 0.04102438815,
                                               0.041907583649999997, 0.05565928034999999,
                                               148.8369807, 271.309704,
                                               127.76516325)
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
            T = 298.0, emax = emax, N = 1, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = False)

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

  params = [0.5872926, 1.0, 0.4240177, 
            0.24294844, 0.85410004, 0.68559576,
            0.70046437, 0.33568583, 0.49608411,
            0.90845018, 0.42513523]
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
  plt.plot(x_sample, res_stress, label = "Sim - %3.0f C interpolate" % (T-273.15))
  plt.plot(df['Nominal_strain'], df['True_stress'], label = "Exp - %3.0f C" % (T-273.15))
  plt.plot(x_sample, stress, label = "Exp - %3.0f C interpolate" % (T-273.15))

  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  # plt.savefig("tension-Ti.png")
  plt.show()
  plt.close()
