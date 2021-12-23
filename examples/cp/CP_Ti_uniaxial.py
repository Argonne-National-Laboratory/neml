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
from tqdm import tqdm


#================================================#
def make_Ti_model(T, emax, Ngrains, erate, nthreads):
#================================================#

  res = CP_Ti_Maker(T = T, emax = emax, N = Ngrains, 
            strain_rate = erate, nthreads = nthreads, 
            verbose = False, Taylor = True,
            PTR = True)
  return res
  
#================================================#
def load_file(path, T):
#================================================#
  fnames = glob.glob(path + "*.csv")
  for f in fnames:
    strain_rate = os.path.basename(f).split('_')[0]
    if strain_rate == "1e-2":
      temp = os.path.basename(f).split('_')[1].split('.')[0]
    if strain_rate == "1e-2" and temp == str(int(T)) + "k":
      df = pd.read_csv(f, usecols=[0,1], names=['True_strain', 'True_stress'], header=None)
      return df



if __name__ == "__main__":

  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Huang-2007-MSEA/"
   # sets up x_scale for both experiment and simulation
  Nsample = 200
  emax = 0.2
  x_sample = np.linspace(0.0, emax*0.99, Nsample)
  #  model grains and threads
  erate = 1.0e-2
  Ngrains = 200
  nthreads = 30 
  # sets up temperature list
  Ts = np.array([298.0, 423.0, 523.0, 623.0, 773.0, 873.0, 973.0])
  
  for T in tqdm(Ts):
    df = load_file(path_1, T)
    res = make_Ti_model(T, emax, Ngrains, erate, nthreads)
    
    true_strain = np.log(1+res['strain'])
    true_stress = res['stress'] * (1 + res['strain'])    
    plt.plot(true_strain, true_stress, 'g--', label = 'fit')
    plt.plot(df['True_strain'], df['True_stress'], 'b-', label = 'data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("CP_Ti_{}".format(int(T)))
    plt.legend()
    plt.grid(True)
    plt.savefig("CP_Ti_calibration_at_{}.png".format(int(T)))
    # plt.show()
    plt.close()
   
 
 
  
