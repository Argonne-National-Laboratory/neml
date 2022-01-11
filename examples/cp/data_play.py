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
T = 1073.0


#================================================#
def load_file(path):
#================================================#
  fnames = glob.glob(path + "*.csv")
  for f in fnames:
    strain_rate = os.path.basename(f).split('_')[0]
    if strain_rate == "1e-2":
      temp = os.path.basename(f).split('_')[1].split('.')[0]
    if strain_rate == "1e-2" and temp == "1073k":
      df = pd.read_csv(f, usecols=[0,1], names=['True_strain', 'True_stress'], header=None)
      return df

if __name__ == "__main__":
 
 
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Yapici-2014-MD/"
  df = load_file(path_1)
  print(df)
  max_stress = df['True_stress'].idxmax()
  new_df = df[:max_stress]

  print(new_df)

















