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
        tau_D1_min, tau_D2_min, tau_D3_min) = (1.6588800000000004, 13.415912979999996,
                                               5.269677479999999, 860.04205774848,
                                               0.09999999999999998, 0.006057790179999999,
                                               0.004386340539999999, 0.335489047695,
                                               11.269026379999998, 25.236950999999994,
                                               814.5508473000001)

max_theta = (X_s_max, k1_1_max, k1_2_max, k1_3_max,
        X_max, g_1_max, g_2_max, g_3_max,
        tau_D1_max, tau_D2_max, tau_D3_max) = (2.0736000000000003, 120.74321681999999,
                                               47.42709732, 1075.0525721856,
                                               0.9, 0.05452011162,
                                               0.03947706486, 0.5032335715425,
                                               101.42123742, 227.132559,
                                               1221.8262709500002)


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
def interpolate_obs(x_sample):
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1)
  stress = interpolate(df['Nominal_strain'], df['True_stress'], x_sample)  
  return stress


if __name__ == "__main__":
  
  T = 298.0

  params = (0.58213844, 0.31971182, 0.0024615, 0.57756845,
            0.54379664, 0.75479085, 0.81933563, 0.81199687,
            0.09674807, 0.41019141, 0.77107571)
          
  (X_s_i, k1_1_i, k1_2_i, k1_3_i, 
        X_i, g_1_i, g_2_i, g_3_i, 
        tau_D1_i, tau_D2_i, tau_D3_i) = convert_to_real(params)

  print(convert_to_real(params))