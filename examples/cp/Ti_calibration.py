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
(taus_1_i, taus_2_i, taus_3_i,
        taut_1_i, taut_2_i, X_s_i,
        k1_1_i, k1_2_i, k1_3_i, 
        X_i, g_1_i, g_2_i, g_3_i, 
        tau_D1_i, tau_D2_i, tau_D3_i) = (170.0, 90.5, 210.0,
                                         180.0, 250.0, 0.9,
                                         1.50, 0.25, 5.00,
                                         0.35, 0.002, 0.002, 0.0055,
                                         100.0, 100.0, 100.0)

sf = 0.9

# sets up x_scale for both experiment and simulation
emax = 0.15
Nsample = 1000
x_sample = np.linspace(0.0, emax*0.99, Nsample)

# sets up the parameters range 
min_theta = (taus_1_min, taus_2_min, taus_3_min, taut_1_min, taut_2_min, 
        X_s_min, k1_1_min, k1_2_min, k1_3_min, 
        X_min, g_1_min, g_2_min, g_3_min,
        tau_D1_min, tau_D2_min, tau_D3_min) = (taus_1_i*(1-sf), taus_2_i*(1-sf), taus_3_i*(1-sf),
                                               taut_1_i*(1-sf), taut_2_i*(1-sf), X_s_i*(1-sf), 
                                               k1_1_i*(1-sf), k1_2_i*(1-sf), k1_3_i*(1-sf), 
                                               X_i*(1-sf), g_1_i*(1-sf),
                                               g_2_i*(1-sf), g_3_i*(1-sf), 
                                               tau_D1_i*(1-sf), tau_D2_i*(1-sf),
                                               tau_D3_i*(1-sf))
                                               
max_theta = (taus_1_max, taus_2_max, taus_3_max, taut_1_max, taut_2_max,
        X_s_max, k1_1_max, k1_2_max, k1_3_max, 
        X_max, g_1_max, g_2_max, g_3_max,
        tau_D1_max, tau_D2_max, tau_D3_max) = (taus_1_i*(1+sf), taus_2_i*(1+sf), taus_3_i*(1+sf),
                                               taut_1_i*(1+sf), taut_2_i*(1+sf), X_s_i*(1+sf), 
                                               k1_1_i*(1+sf), k1_2_i*(1+sf), k1_3_i*(1+sf), 
                                               X_i*(1+sf), g_1_i*(1+sf),
                                               g_2_i*(1+sf), g_3_i*(1+sf), 
                                               tau_D1_i*(1+sf), tau_D2_i*(1+sf),
                                               tau_D3_i*(1+sf))

#================================================#
def convert_to_real(p):
#================================================#
    bounds = np.array([
                    [taus_1_min, taus_1_max],    # tau_s1
                    [taus_2_min, taus_2_max],    # tau_s2
                    [taus_3_min, taus_3_max],    # tau_s3
                    [taut_1_min, taut_1_max],    # tau_t1
                    [taut_2_min, taut_2_max],    # tau_t2
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
         params[10], params[11], params[12],
         params[13], params[14], params[15])
  
  theta = convert_to_real(theta_in)
  
  taus_1, taus_2, taus_3, taut_1, taut_2, X_s, k1_1, k1_2, k1_3, X, g_1, g_2, g_3, tau_D1, tau_D2, tau_D3  = (theta[0],
        theta[1],theta[2],theta[3],theta[4],theta[5],theta[6],theta[7],theta[8],theta[9],theta[10],
        theta[11],theta[12],theta[13],theta[14],theta[15])
  
  res = make_model(taus_1, taus_2, taus_3, taut_1, taut_2,
            X_s, k1_1, k1_2, k1_3, X, 
            g_1, g_2, g_3,
            tau_D1, tau_D2, tau_D3,
            T = 298.0, emax = emax, N = 50, 
            strain_rate = 1.0e-4, nthreads = 30, 
            verbose = False, Taylor = True,
            PTR = True)
  return res


#================================================#
def interpolate_obs():
#================================================#
  # interpolate real experimental data
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"
  df = load_file(path_1) 
  return df

#================================================#
def R(params):
#================================================#
  res = make_Ti_model(params)
  df = interpolate_obs()
  sim = interpolate(df['Nominal_strain'], df['True_stress'], x_sample)
  yobs = interpolate(res['strain'], res['stress'], x_sample)
  R = la.norm(sim - yobs)
  print("Current residual: %e" % R)
  return R



#================================================#
def set_scale(p_theta, min_theta, max_theta, control_param = False):
#================================================#
  theta_in  = (p_theta[0],p_theta[1],p_theta[2],p_theta[3],p_theta[4],p_theta[5],p_theta[6],
               p_theta[7], p_theta[8], p_theta[9], p_theta[10], p_theta[11], p_theta[12],
               p_theta[13], p_theta[14], p_theta[15])
  min_in  = np.array([min_theta[0],min_theta[1],min_theta[2],min_theta[3],min_theta[4],
                      min_theta[5],min_theta[6],min_theta[7],min_theta[8],min_theta[9],
                      min_theta[10],min_theta[11],min_theta[12],min_theta[13],
                      min_theta[14],min_theta[15]])
  max_in = np.array([max_theta[0],max_theta[1],max_theta[2],max_theta[3],max_theta[4],
                     max_theta[5],max_theta[6],max_theta[7],max_theta[8],max_theta[9],
                     max_theta[10],max_theta[11],max_theta[12],max_theta[13],
                     max_theta[14],max_theta[15]])

  if control_param:
    for i in range(len(theta_in)):
      if theta_in[i] == 0.0:
        if i == 5 or i == 9:
          min_in[i] = 0.1
          max_in[i] = 1.0
          for j in range(len(theta_in)):
            if j != 5 and j != 9:
              min_in[j] = min_in[j] * 0.8
              max_in[j] = max_in[j] * 1.5            
        else:
          max_in[i] = min_in[i]
          min_in[i] = min_in[i] * 0.8
      elif theta_in[i] == 1.0:
        if i == 5 or i == 9:
          min_in[i] = 0.1
          max_in[i] = 1.0
          for j in range(len(theta_in)):
            if j != 5 and j != 9:
              min_in[j] = min_in[j] * 0.8
              max_in[j] = max_in[j] * 1.5            
        else:
          min_in[i] = max_in[i]
          max_in[i] = max_in[i] * 1.5
      else:
        print('number:', i, 'is fine!!')
  else:
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
  
  # select optimize mode
  easy_mode = True
  
  # here we go
  taus_1_range = [0.0, 1.0]
  taus_2_range = [0.0, 1.0]
  taus_3_range = [0.0, 1.0]
  taut_1_range = [0.0, 1.0]
  taut_2_range = [0.0, 1.0]
  X_s_range = [0.0, 1.0]
  k1_1_range = [0.0, 1.0]
  k1_2_range = [0.0, 1.0]
  k1_3_range = [0.0, 1.0]
  X_range = [0.0, 1.0]
  g_1_range = [0.0, 1.0]
  g_2_range = [0.0, 1.0]
  g_3_range = [0.0, 1.0]
  tau_D1_range = [0.0, 1.0]
  tau_D2_range = [0.0, 1.0]
  tau_D3_range = [0.0, 1.0]
  
  
  p0 = [ra.uniform(*taus_1_range), ra.uniform(*taus_2_range), ra.uniform(*taus_3_range), 
       ra.uniform(*taut_1_range), ra.uniform(*taut_2_range),
       ra.uniform(*X_s_range), ra.uniform(*k1_1_range),
       ra.uniform(*k1_2_range), ra.uniform(*k1_3_range),
       ra.uniform(*X_range), ra.uniform(*g_1_range),
       ra.uniform(*g_2_range), ra.uniform(*g_3_range),
       ra.uniform(*tau_D1_range), ra.uniform(*tau_D2_range),
       ra.uniform(*tau_D3_range)]
  

  if easy_mode:
    res = opt.minimize(R, p0, method = 'L-BFGS-B', bounds = [taus_1_range, taus_2_range, taus_3_range,
        taut_1_range, taut_2_range, X_s_range,
        k1_1_range, k1_2_range, k1_3_range,
        X_range, g_1_range, g_2_range,
        g_3_range, tau_D1_range,
        tau_D2_range, tau_D3_range])

    print(res.success)
    print(res)

    ref_n = (taus_1_n, taus_2_n, taus_3_n, taut_1_n, taut_2_n, 
             X_s_n, k1_1_n, k1_2_n, k1_3_n, X_n, g_1_n, g_2_n, g_3_n, 
             tau_D1_n, tau_D2_n, tau_D3_n) = (res.x[0], res.x[1], res.x[2], 
                                              res.x[3], res.x[4], res.x[5], 
                                              res.x[6], res.x[7], res.x[8],
                                              res.x[9], res.x[10], res.x[11],
                                              res.x[12], res.x[13], res.x[14],
                                              res.x[15])

    res_final = make_Ti_model(ref_n)
    df = interpolate_obs()
    # visualize the fitting
    T = 298.0
    plt.plot(res_final['strain'], res_final['stress'], label = "Sim - %3.0f C" % (T-273.15))
    plt.plot(df['Nominal_strain'], df['True_stress'], label = "Exp - %3.0f C" % (T-273.15))
    plt.legend(loc='best')
    plt.xlabel("Strain (mm/mm)")
    plt.ylabel("Stress (MPa)")
    plt.savefig("tension-Ti.png")
    plt.show()
    plt.close()

  else:
  
    flag = False
    iq = 0
    while not flag:
      res = minimize_parallel(R, p0, bounds = [taus_1_range, taus_2_range, taus_3_range,
        taut_1_range, taut_2_range, X_s_range,
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

    ref_n = (taus_1_n, taus_2_n, taus_3_n, taut_1_n, taut_2_n, 
             X_s_n, k1_1_n, k1_2_n, k1_3_n, X_n, g_1_n, g_2_n, g_3_n, 
             tau_D1_n, tau_D2_n, tau_D3_n) = (res.x[0], res.x[1], res.x[2], 
                                              res.x[3], res.x[4], res.x[5], 
                                              res.x[6], res.x[7], res.x[8],
                                              res.x[9], res.x[10])

    switch = False
 
    while not switch:
      dex = 0
      for i in range(len(ref_n)):
        if ref_n[i] <= 0.0 or ref_n[i] >= 1.0:
          switch = False
          flag = False
          min_out, max_out = set_scale(ref_n, min_theta, max_theta, control_param = False)
          min_theta = (taus_1_min, taus_2_min, taus_3_min, taut_1_min, taut_2_min, 
                       X_s_min, k1_1_min, k1_2_min, k1_3_min, 
                       X_min, g_1_min, g_2_min, g_3_min,
                       tau_D1_min, tau_D2_min, tau_D3_min) = (min_out[0],min_out[1],
                                                              min_out[2],min_out[3],
                                                              min_out[4],min_out[5],
                                                              min_out[6],min_out[7],
                                                              min_out[8],min_out[9],
                                                              min_out[10],min_out[11],
                                                              min_out[12],min_out[13],
                                                              min_out[14],min_out[15])
                    
          max_theta = (taus_1_max, taus_2_max, taus_3_max, taut_1_max, taut_2_max,
                       X_s_max, k1_1_max, k1_2_max, k1_3_max, 
                       X_max, g_1_max, g_2_max, g_3_max,
                       tau_D1_max, tau_D2_max, tau_D3_max) = (max_out[0],max_out[1],
                                                              max_out[2],max_out[3],
                                                              max_out[4],max_out[5],
                                                              max_out[6],max_out[7],
                                                              max_out[8],max_out[9],
                                                              max_out[10],max_out[11],
                                                              max_out[12],max_out[13],
                                                              max_out[14],max_out[15])
                
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
        print('min=:', taus_1_min, taus_2_min, taus_3_min, taut_1_min, taut_2_min, 
                       X_s_min, k1_1_min, k1_2_min, k1_3_min, 
                       X_min, g_1_min, g_2_min, g_3_min,
                       tau_D1_min, tau_D2_min, tau_D3_min)
                       
        print('max=:', taus_1_max, taus_2_max, taus_3_max, taut_1_max, taut_2_max,
                       X_s_max, k1_1_max, k1_2_max, k1_3_max, 
                       X_max, g_1_max, g_2_max, g_3_max,
                       tau_D1_max, tau_D2_max, tau_D3_max)
            
      while not flag:
        res = minimize_parallel(R, p0, bounds = [taus_1_range, taus_2_range, taus_3_range,
                                                 taut_1_range, taut_2_range, X_s_range,
                                                 k1_1_range, k1_2_range, k1_3_range,
                                                 X_range, g_1_range, g_2_range,
                                                 g_3_range, tau_D1_range,
                                                 tau_D2_range, tau_D3_range])
        print(res.success)
        if res.success == True:
            flag = True
            print(res.success)
            print(res)
            ref_n = (taus_1_n, taus_2_n, taus_3_n, taut_1_n, taut_2_n, 
                     X_s_n, k1_1_n, k1_2_n, k1_3_n, X_n, g_1_n, g_2_n, g_3_n, 
                     tau_D1_n, tau_D2_n, tau_D3_n) = (res.x[0], res.x[1], res.x[2], 
                                                      res.x[3], res.x[4], res.x[5],
                                                      res.x[6], res.x[7], res.x[8],
                                                      res.x[9], res.x[10], res.x[11],
                                                      res.x[12], res.x[13], res.x[14],
                                                      res.x[15])
            dex = 0
        iq += 1
        if iq >=50:
          raise ValueError("Not able to optimize the initialize > 50") 
          
    res_final = make_Ti_model(ref_n)
    df = interpolate_obs()
    # visualize the fitting
    T = 298.0
    plt.plot(res_final['strain'], res_final['stress'], label = "Sim - %3.0f C" % (T-273.15))
    plt.plot(df['Nominal_strain'], df['True_stress'], label = "Exp - %3.0f C" % (T-273.15))
    plt.legend(loc='best')
    plt.xlabel("Strain (mm/mm)")
    plt.ylabel("Stress (MPa)")
    # plt.savefig("tension-Ti.png")
    plt.show()
    plt.close()