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


# sets up initial model parameters
(X_s_i, k1_1_i, k1_2_i, k1_3_i, 
        X_i, g_1_i, g_2_i, g_3_i, 
        tau_D1_i, tau_D2_i, tau_D3_i) = (0.9, 75.0, 12.5, 250.0,
                                         0.3, 0.02, 0.02, 0.055,
                                         100.0, 100.0, 100.0)
                                         
sf = 0.2

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

  res_stress = interpolate(res['strain'], res['stress'], x_sample) 
  return res_stress


#================================================#
def interpolate_obs(x_sample):
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1)
  stress = interpolate(df['Nominal_strain'], df['True_stress'], x_sample)  
  # stress = inter.interp1d(df['Nominal_strain'], df['True_stress'])(res['strain']) 
  return stress

#================================================#
def R(params):
#================================================#

  res = make_Ti_model(params)
  yobs = interpolate_obs(x_sample)
  R = la.norm(yobs - res)
  print("Current residual: %e" % R)
  return R



#================================================#
def set_scale(p_theta, min_theta, max_theta):
#================================================#
  theta_in  = (p_theta[0],p_theta[1],p_theta[2],p_theta[3],p_theta[4],p_theta[5],p_theta[6],
               p_theta[7], p_theta[8], p_theta[9], p_theta[10])
  min_in  = np.array([min_theta[0],min_theta[1],min_theta[2],min_theta[3],min_theta[4],
                      min_theta[5],min_theta[6],min_theta[7],min_theta[8],min_theta[9],
                      min_theta[10]])
  max_in = np.array([max_theta[0],max_theta[1],max_theta[2],max_theta[3],max_theta[4],
                     max_theta[5],max_theta[6],max_theta[7],max_theta[8],max_theta[9],
                     max_theta[10]])
  for i in range(len(theta_in)):
    if theta_in[i] == 0.0:
      max_in[i] = min_in[i]
      min_in[i] = min_in[i] * 0.8
    elif theta_in[i] == 1.0:
      min_in[i] = max_in[i]
      max_in[i] = max_in[i] * 1.5
    else:
      print('number:', i, 'is fine!!')
            
  min_out  = min_in
  max_out = max_in
    
  return min_out, max_out


if __name__ == "__main__":
  
  # here we go
  rranges = (slice(0, 1, 0.01), slice(0, 1, 0.01), slice(0, 1, 0.01),
                slice(0, 1, 0.01), slice(0, 1, 0.01), slice(0, 1, 0.01),
                slice(0, 1, 0.01), slice(0, 1, 0.01), slice(0, 1, 0.01),
                slice(0, 1, 0.01))

  flag = False
  iq = 0
  while not flag:
    res = opt.brute(R, rranges, full_output=True,
                          finish=opt.fmin)
    print(res.success)
    if res.success == True:
      flag = True
      print(res.success)
      print(res.x)
    iq += 1
    if iq >=3:
      raise ValueError("Not able to optimize the initialize > 3")
  
