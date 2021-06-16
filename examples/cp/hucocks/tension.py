#!/usr/bin/env python3

from generate_model import *

from neml import drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 50
  nthreads = 3
  model = make_model(N, nthreads = nthreads)
 
  Ts = np.array([500,550,600,650.0]) + 273.15

  for T in Ts:
    res = drivers.uniaxial_test(model, 8.33e-5, T = T, verbose = True)
  
    plt.plot(res['strain'], res['stress'], label = "%3.0f C" % (T-273.15))

  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.savefig("tension-unaged.png")

