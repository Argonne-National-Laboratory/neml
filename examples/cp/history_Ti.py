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

Ngrains = 50
nthreads = 1



#================================================#
def Ti_ld_maker(taus_1, taus_2, taus_3,
                taut_1, taut_2, X_s,
                k1_1, k1_2, k1_3,
                k2_1, k2_2, k2_3,
                N = Ngrains, nthreads = nthreads,
                large_deform = False):
#================================================#  

  res = Ti_maker_Polycrystal(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s,
            k1_1, k1_2, k1_3,
            k2_1, k2_2, k2_3,
            T = 296.0, emax = 0.05, N = 1, steps = 100,
            strain_rate = 1.0e-4, nthreads = 1,
            verbose = True, Taylor = True,
            PTR = True, return_hardening = False,
            full_results = False,
            large_deform = large_deform)
  return res


#================================================#
def make_simple_cubic(N = Ngrains, nthreads = nthreads):
#================================================#  
  smodel = tensile_single_cubic(simplelinear = True,         
            verbose = False,
            update_rotation = True)

  orientations = rotations.random_orientations(N)

  model = polycrystal.TaylorModel(smodel, orientations, nthreads = nthreads)

  return model

#================================================#
def make_simple_Ti(taus_1, taus_2, taus_3,
            taut_1, taut_2, H1, H2, 
            N = Ngrains, nthreads = nthreads):
#================================================#  
  smodel = make_simple_singlecrystal(taus_1, taus_2, taus_3,
            taut_1, taut_2, H1, H2,
            verbose = True, PTR = False, 
            return_hardening = False,
            update_rotation = False)

  orientations = rotations.random_orientations(N)

  model = polycrystal.TaylorModel(smodel, orientations, nthreads = nthreads)

  return model
#================================================#
def make_Ti_polycrystal(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s,
            k1_1, k1_2, k1_3,
            k2_1, k2_2, k2_3, 
            N = Ngrains, nthreads = nthreads):
#================================================#  
  smodel = make_Ti_singlecrystal(taus_1, taus_2, taus_3,
            taut_1, taut_2, X_s,
            k1_1, k1_2, k1_3,
            k2_1, k2_2, k2_3,
            verbose = True, PTR = False, 
            return_hardening = False,
            update_rotation = False)

  orientations = rotations.random_orientations(N)

  model = polycrystal.TaylorModel(smodel, orientations, nthreads = nthreads)

  return model
#================================================#
def load_file(path):
#================================================#
  fnames = glob.glob(path + "*.csv")
  for f in fnames:
    strain_rate = os.path.basename(f).split('_')[0]
    if strain_rate == "1e-2":
      temp = os.path.basename(f).split('_')[1].split('.')[0]
    if strain_rate == "1e-2" and temp == "298k":
      df = pd.read_csv(f, usecols=[0,1], names=['True_strain', 'True_stress'], header=None)
      return df

if __name__ == "__main__":
 
  # define the model trying to use
  use_model = "Ti_polycrystal"
  # sets up x_scale for both experiment and simulation
  emax = 0.2
  Nsample = 200
  x_sample = np.linspace(0.0, emax*0.99, Nsample)
  #  model grains and threads
  erate = 1.0e-2
  T = 298.0
 
  path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Huang-2007-MSEA/"
  df = load_file(path_1)
  print(df)
  
  
  if use_model == "Ti_polycrystal":
    # room temperature
    params = [190.0, 110.5, 230.0,
        180.0, 250.0, 0.9,
        1.0, 0.25, 5.0,
        25.0, 25.0, 25.0]
    tmodel = make_Ti_polycrystal(*params)
  elif use_model == "Ti_simple":
    # 973k temperature
    params = [25.0, 25.0, 25.0,
        100.0, 180.0, 1.0, 1.0]
    tmodel = make_simple_Ti(*params)
  elif use_model == "cubic":
    tmodel = make_simple_cubic()
  
  
  # output histroy internal variables
  full_results = False
  
  if full_results:
    res = drivers.uniaxial_test(tmodel, erate = erate, emax = emax,
                              T = T, verbose = True, full_results = True)

    density = np.array(res['history'])
    length = 47 * Ngrains # every single grain will have size of 47 for the internal history vector
    starting = np.arange(8, length, 32)  # every 32 colomns will store the history variables
    # store the dislocation density evolution for each grains
    dd = np.vstack((density[:,i:i+12].T for i in starting[:-1])).T
    df = pd.DataFrame(density)
    df.to_csv('histroy_checkup_{}.csv'.format(int(T)))
    dd_df = pd.DataFrame(dd)
    dd_df.to_csv('dd_checkup_{}.csv'.format(int(T)))

    true_strain = np.log(1+np.abs(res['strain']))
    true_stress = np.array(res['stress']) * (1 + np.abs(res['strain']))

    print("The dislocation evolution is:")
    forest_density = np.mean(dd, axis=1) * 1.0e18
    true_strain = true_strain[:, 2]
    true_stress = np.abs(true_stress[:, 2])


    plt.plot(true_strain, forest_density, 'g--', label = 'forest dislocation density')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    # plt.savefig("Dislocation-history-{}.png".format(int(T)))
    plt.show()
    plt.close()
  else:
    res = drivers.uniaxial_test(tmodel, erate = erate, emax = emax,
                              T = T, verbose = True, full_results = full_results)
       
    true_strain = np.log(1+np.abs(res['strain']))
    true_stress = np.array(res['stress']) * (1 + np.abs(res['strain']))     
    plt.plot(true_strain, true_stress, 'g--', label = 'stress-strain')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    # plt.savefig("Dislocation-history-{}.png".format(int(T)))
    plt.show()
    plt.close()    
  
  
