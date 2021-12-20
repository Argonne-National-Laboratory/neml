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
erate = 1.0e-4
Ngrains = 10
nthreads = 1
T = 296.0

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
            verbose = True, Taylor = True,
            PTR = True)

  return res

#================================================#
def load_file(path, temper):
#================================================#
  for file in glob.glob(path + temper + "k.csv"):
    df = pd.read_csv(file, usecols=[0,1], names=['Nominal_strain', 'True_stress'], header=None)
  return df

  
if __name__ == "__main__":
 
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Nasser-1999-Acta/"
  df = load_file(path_1, str(int(T)))
  
  params = [100.0, 80.5, 100.0,
        130.0, 200.0, 0.9,
        1.00, 0.15, 4.00,
        16.0, 16.0, 16.0]
  res = make_Ti_model(x_sample, *params)
  plt.plot(res['strain'], res['stress'], 'g--', label = 'fit')
  plt.plot(df['Nominal_strain'], df['True_stress'], 'b-', label = 'data')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.legend()
  # plt.savefig("Ti-fitting-simple.png")
  plt.show()
  plt.close()
    

  
