#!/usr/bin/env python3

from generate_model import *

from neml import drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 50
  nthreads = 3
  model = make_model(N, nthreads = nthreads)
 
  Ts = np.array([500,550,600,650.0]) + 273.15
  time = 10000.0 * 3600.0
  lrate = 1.0
  stresses = np.array([40.0,40.0,40.0,40.0])

  for T,s in zip(Ts, stresses):
    print(T)

    res = drivers.creep(model, s, lrate, time, T = T, nsteps = 200, 
        nsteps_up = 20, logspace = True, verbose = True)

    plt.semilogx(res['time']/3600.0, res['strain'], 
        label = "%3.0f C" % (T-273.15))
  
  plt.legend(loc="best")
  plt.xlabel("Time (hr)")
  plt.ylabel("Strain (mm/mm)")
  plt.savefig("creep.png")

