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
def load_file(path, T, rate):
#================================================#
  fnames = glob.glob(path + "*.csv")
  for f in fnames:
    strain_rate = os.path.basename(f).split('_')[0]
    if strain_rate == rate:
      temp = os.path.basename(f).split('_')[1].split('.')[0]
    if strain_rate == rate and temp == str(int(T)) + 'k':
      df = pd.read_csv(f, usecols=[0,1], names=['True_strain', 'True_stress'], header=None)
      return df

if __name__ == "__main__":

  # set up model grains and threads
  Ngrains = 50
  nthreads = 10
  # tensile conditions
  rate = "1e-2"
  emax = 0.2
  erate = float(rate)
  Ts = np.array([298.0, 423.0, 523.0, 623.0, 773.0, 873.0, 973.0, 1073.0, 1173.0])
  
  path_1 = xxxxx
  path_2 = xxxxx

  for T in Ts:
    if T < 1000.0:
      path = path_1
    else:
      path = path_2
    df = load_file(path, T, rate)
    tmodel = make_Ti_polycrystal(Ngrains, nthreads)
    res = drivers.uniaxial_test(tmodel, erate = erate, emax = emax,
                          sdir = np.array([0,0,-1,0,0,0]),
                          T = T, verbose = True, full_results = False)
    # save the stress-strain data for future use
    data = pd.DataFrame({
      "strain": res['strain'],
      "stress": res['stress']}) 
    data.to_csv('res_{}.csv'.format(int(T)))
    # convert to engineering stress-strain
    eng_strain = np.exp(df['True_strain']) - 1.0
    eng_stress = df['True_stress'] / (1 + eng_strain)
    plt.plot(eng_strain, eng_stress, label = 'Exp-{}'.format(int(T)))
    plt.plot(res['strain'], res['stress'], label = 'Model-{}'.format(int(T)))
    plt.xlabel('Strain (mm/mm)')
    plt.ylabel('Stress (MPa)')
    plt.legend()
    plt.grid(True)
    plt.savefig("Engineering-stress-strain-{}.png".format(int(T)))
    plt.show()
    plt.close()    
  
  
