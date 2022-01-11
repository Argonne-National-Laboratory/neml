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
emax = 0.2
Nsample = 200
x_sample = np.linspace(0.0, emax*0.99, Nsample)
#  model grains and threads
erate = 1.0e-2
Ngrains = 500
nthreads = 1
T = 973.0

#================================================#
def make_Ti_model(x_sample, taus_1, taus_2, taus_3,
                  taut_1, taut_2, X_s,
                  k1_1, k1_2, k1_3,
                  k2_1, k2_2, k2_3):
#================================================#

  res = Ti_maker_sim(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s,
            k1_1, k1_2, k1_3,
            k2_1, k2_2, k2_3,
            T = T, emax = emax, N = Ngrains,
            strain_rate = erate, nthreads = nthreads,
            verbose = True, Taylor = False,
            PTR = True, return_hardening = True,
            full_results = True)

  return res
  
#================================================#
def load_file(path):
#================================================#
  fnames = glob.glob(path + "*.csv")
  for f in fnames:
    strain_rate = os.path.basename(f).split('_')[0]
    if strain_rate == "1e-2":
      temp = os.path.basename(f).split('_')[1].split('.')[0]
    if strain_rate == "1e-2" and temp == "973k":
      df = pd.read_csv(f, usecols=[0,1], names=['True_strain', 'True_stress'], header=None)
      return df



if __name__ == "__main__":
 
 
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Huang-2007-MSEA/"
  df = load_file(path_1)
  print(df)
  
  # room temperature
  # 973k temperature
  params = [25.0, 25.0, 25.0,
        100.0, 180.0, 0.9,
        1.0, 0.25, 5.0,
        1000.0, 1000.0, 1000.0]
        
  res, tau_model = make_Ti_model(x_sample, *params)

  
  density = np.array(res['history'])
  
  
  true_strain = np.log(1+np.abs(res['strain']))
  true_stress = np.array(res['stress']) * (1 + np.abs(res['strain']))
  
  print("density shape: ", np.shape(density))
  print("strain shape: ", np.shape(true_strain))
  
  
  density = np.mean(density[:, -12:], axis=1) * 1.0e18
  true_strain = true_strain[:, 2]
  """
  plt.plot(true_strain, true_stress, 'g--', label = 'fit')
  plt.plot(df['True_strain'], df['True_stress'], 'b-', label = 'data')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.legend()
  plt.grid(True)
  # plt.savefig("Ti-fitting-simple.png")
  plt.show()
  plt.close()  
  """
  
  plt.plot(true_strain, density, 'g--', label = 'forest dislocation density')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.legend()
  plt.grid(True)
  plt.savefig("Dislocation-history.png")
  plt.show()
  plt.close()  
  
  
