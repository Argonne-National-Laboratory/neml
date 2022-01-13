#!/usr/bin/env python3

from Ti_maker import *
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
def make_Ti_polycrystal(N, nthreads):
#================================================#  
  smodel = Ti_singlecrystal(verbose = True, PTR = True, 
                            return_hardening = False,
                            update_rotation = True)

  orientations = rotations.random_orientations(N)

  model = polycrystal.TaylorModel(smodel, orientations, nthreads = nthreads)

  return model
#================================================#
def load_file(path, T):
#================================================#
  fnames = glob.glob(path + "*.csv")
  for f in fnames:
    strain_rate = os.path.basename(f).split('_')[0]
    if strain_rate == "1e-2":
      temp = os.path.basename(f).split('_')[1].split('.')[0]
    if strain_rate == "1e-2" and temp == str(int(T)) + 'k':
      df = pd.read_csv(f, usecols=[0,1], names=['True_strain', 'True_stress'], header=None)
      return df

if __name__ == "__main__":

  Ngrains = 50
  nthreads = 1
  # sets up x_scale for both experiment and simulation
  emax = 0.2
  erate = 1.0e-2
  Ts = np.array([298.0, 423.0, 523.0, 623.0, 773.0, 873.0, 973.0, 1073.0, 1173.0])
  
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Huang-2007-MSEA/"
  path_2 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Yapici-2014-MD/"
  
  for T in Ts:
    if T < 1000.0:
      path = path_1
    else:
      path = path_2
    df = load_file(path, T)
    tmodel = make_Ti_polycrystal(Ngrains, nthreads)

    res = drivers.uniaxial_test(tmodel, erate = erate, emax = emax,
                          T = T, verbose = True, full_results = False)

    eng_strain = np.exp(df['True_strain']) - 1.0
    eng_stress = df['True_stress'] / (1 + eng_strain)
    plt.plot(eng_strain, eng_stress, label = 'Exp')
    plt.plot(res['strain'], res['stress'], label = 'Ti Model')
    plt.xlabel('Strain (mm/mm)')
    plt.ylabel('Stress (MPa)')
    plt.legend()
    plt.grid(True)
    plt.savefig("Engineering-stress-strain-{}.png".format(int(T)))
    plt.show()
    plt.close()    
  
  
