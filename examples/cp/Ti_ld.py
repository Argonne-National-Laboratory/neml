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

from neml.cp import hucocks, crystallography, sliprules, slipharden, inelasticity, kinematics, singlecrystal, polycrystal
from neml.math import rotations


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
 
  # sets up x_scale for both experiment and simulation
  Ngrains = 5
  Nthreads = 10
  emax = 0.2
  Nsample = 200
  x_sample = np.linspace(0.0, emax*0.99, Nsample)
  #  model grains and threads
  erate = 1.0e-2
  T = 973.0
  steps = 100
 
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Huang-2007-MSEA/"
  df = load_file(path_1)

  # room temperature
  params = [190.0, 110.5, 230.0,
        180.0, 250.0, 0.9,
        1.0, 0.25, 5.0,
        250.0, 250.0, 250.0]

  res = Ti_maker_Polycrystal(*params,
            T = T, emax = emax, N = Ngrains, steps = steps,
            strain_rate = erate, nthreads = Nthreads,
            verbose = True, Taylor = True,
            PTR = True, return_hardening = False,
            full_results = False,
            large_deform = True)
  
  plt.plot(res[0], res[1], label = 'Large deformation')
  plt.xlabel('Strain (mm/mm)')
  plt.ylabel('Stress (MPa)')
  plt.legend()
  plt.grid(True)
  # plt.savefig("Dislocation-history-{}.png".format(int(T)))
  plt.show()
  plt.close()