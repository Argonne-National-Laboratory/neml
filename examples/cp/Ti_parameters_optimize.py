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



# sets up the parameters range 
min_theta = (X_s_min, k1_1_min, k1_2_min, k1_3_min, 
        X_min, g_1_min, g_2_min, g_3_min,
        tau_D1_min, tau_D2_min, tau_D3_min) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                               0.0, 0.0, 0.0, 0.0, 0.0)
                                               
max_theta = (X_s_max, k1_1_max, k1_2_max, k1_3_max, 
        X_max, g_1_max, g_2_max, g_3_max,
        tau_D1_max, tau_D2_max, tau_D3_max) = (1.0, 1.5, 2.5, 5.0, 0.5,
                                                0.1, 0.1, 0.1, 200.0, 200.0,
                                                200.0)


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
            T = 298.0, N = 1, 
            strain_rate = 1.0e-4, nthreads = 1, 
            verbose = True)

  return res


#================================================#
def interpolate_obs(res):
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1)
  # stress = interpolate(df['Nominal_strain'], df['True_stress'], res['strain'])  
  print(res['strain'])
  stress = inter.interp1d(df['Nominal_strain'], df['True_stress'])(res['strain']) 
  return stress

#================================================#
def R(params):
#================================================#

  res = make_Ti_model(params)
  yobs = interpolate_obs(res)
  R = la.norm(yobs - res['stress'])
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
  X_s_range = [1.0e-3, 1.0]
  k1_1_range = [1.0e-3, 1.0]
  k1_2_range = [1.0e-3, 1.0]
  k1_3_range = [1.0e-3, 1.0]
  X_range = [1.0e-3, 1.0]
  g_1_range = [1.0e-3, 1.0]
  g_2_range = [1.0e-3, 1.0]
  g_3_range = [1.0e-3, 1.0]
  tau_D1_range = [1.0e-3, 1.0]
  tau_D2_range = [1.0e-3, 1.0]
  tau_D3_range = [1.0e-3, 1.0]
  
  
  p0 = [ra.uniform(*X_s_range), ra.uniform(*k1_1_range),
       ra.uniform(*k1_2_range), ra.uniform(*k1_3_range),
       ra.uniform(*X_range), ra.uniform(*g_1_range),
       ra.uniform(*g_2_range), ra.uniform(*g_3_range),
       ra.uniform(*tau_D1_range), ra.uniform(*tau_D2_range),
       ra.uniform(*tau_D3_range)]
  

  flag = False
  iq = 0
  while not flag:
    res = minimize_parallel(R, p0, bounds = [X_s_range,
        k1_1_range, k1_2_range, k1_3_range,
        X_range, g_1_range, g_2_range,
        g_3_range, tau_D1_range,
        tau_D2_range, tau_D3_range])
    print(res.success)
    if res.success == True:
      flag = True
      print(res.success)
      print(res.x)
    iq += 1
    if iq >=3:
      raise ValueError("Not able to optimize the initialize > 3")

  ref_n = (X_s_n, k1_1_n, k1_2_n, k1_3_n, X_n, g_1_n, g_2_n, g_3_n, tau_D1_n, tau_D2_n, tau_D3_n) = (res.x[0], 
                res.x[1], res.x[2], res.x[3], res.x[4],
                res.x[5], res.x[6], res.x[7], res.x[8],
                res.x[9], res.x[10])
  
  switch = False
 
  while not switch:
    dex = 0
    for i in range(len(ref_n)):
      if ref_n[i] <= 0.0 or ref_n[i] >= 1.0:
        switch = False
        flag = False
        min_out, max_out = set_scale(ref_n, min_theta, max_theta)
        min_theta = (X_s_min, k1_1_min, k1_2_min, k1_3_min, 
            X_min, g_1_min, g_2_min, g_3_min,
            tau_D1_min, tau_D2_min, tau_D3_min) = (min_out[0],min_out[1],
                    min_out[2],min_out[3],min_out[4],min_out[5],min_out[6],
                    min_out[7],min_out[8],min_out[9],min_out[10])
                    
        max_theta = (X_s_max, k1_1_max, k1_2_max, k1_3_max, 
            X_max, g_1_max, g_2_max, g_3_max,
            tau_D1_max, tau_D2_max, tau_D3_max) = (max_out[0],max_out[1],
                    max_out[2],max_out[3],max_out[4],max_out[5],max_out[6],
                    max_out[7],max_out[8],max_out[9],max_out[10])
                
        dex += 0
      else:
        dex += 1
            
    if dex < len(ref_n) :
      switch = False
      flag = False
    else:
      switch = True
      flag = True
      print('all the params are fine now!!')
      print('min=:', X_s_min, k1_1_min, k1_2_min, k1_3_min, 
            X_min, g_1_min, g_2_min, g_3_min,
            tau_D1_min, tau_D2_min, tau_D3_min)
      print('max=:', X_s_max, k1_1_max, k1_2_max, k1_3_max, 
            X_max, g_1_max, g_2_max, g_3_max,
            tau_D1_max, tau_D2_max, tau_D3_max)
            
    while not flag:
      res = minimize_parallel(R, p0, bounds = [X_s_range,
        k1_1_range, k1_2_range, k1_3_range,
        X_range, g_1_range, g_2_range,
        g_3_range, tau_D1_range,
        tau_D2_range, tau_D3_range])
      print(res.success)
      if res.success == True:
          flag = True
          print(res.success)
          print(res.x)
          ref_n = (X_s_n, k1_1_n, k1_2_n, k1_3_n, X_n, g_1_n, g_2_n, g_3_n, tau_D1_n, tau_D2_n, tau_D3_n) = (res.x[0], 
              res.x[1], res.x[2], res.x[3], res.x[4],
              res.x[5], res.x[6], res.x[7], res.x[8],
              res.x[9], res.x[10])
          dex = 0
      iq += 1
      if iq >=50:
        raise ValueError("Not able to optimize the initialize > 50") 