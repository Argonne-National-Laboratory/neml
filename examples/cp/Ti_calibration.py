#!/usr/bin/env python3

import sys
import os, glob
import numpy as np
import numpy.random as ra

import scipy.interpolate as inter

from generate_model_Ti import *
from neml import drivers
import matplotlib.pyplot as plt

import pandas as pd
import xarray as xr
import tqdm
import warnings
warnings.filterwarnings("ignore")

def simulation(strain_rate, N, nthreads, T):
  
  model = make_model(strain_rate, N, nthreads = nthreads)
  res = drivers.uniaxial_test(model, strain_rate, T = T, verbose = False)

  return res


def interpolate(strain, stress, targets):
  """
    This is just to make sure all our values line up, in case the model
    adaptively integrated or something
  """
  return inter.interp1d(strain, stress)(targets) 


def load_file(path):
  for file in glob.glob(path + "CG_Ti.csv"):
    df = pd.read_csv(file, usecols=[0,1], names=['Nominal_strain', 'True_stress'], header=None)
  return df

path_1 = "/mnt/c/Users/ladmin/Desktop/argonne/RTRC_data_extract/Ito-2019-MSEB/"


if __name__ == "__main__":
  N = 100
  nthreads = 1
  strain_rate = 1.0e-4
  T = 298.0
  # samples = 200

  res = simulation(strain_rate, N, nthreads, T)
  df = load_file(path_1)
  
  stress = interpolate(df['Nominal_strain'], df['True_stress'], res['strain'])
  
  plt.plot(res['strain'], res['stress'], label = "Sim - %3.0f C" % (T-273.15))
  plt.plot(res['strain'], stress, label = "Exp - %3.0f C" % (T-273.15))

  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  # plt.savefig("tension-Ti.png")
  plt.show()
  plt.close()