#!/usr/bin/env python3

from generate_model import *

from neml import drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 50
  nthreads = 1
  model = make_model(N, nthreads = nthreads)

  aging = 100000.0*3600.0
  nage = 50
 
  Ts = np.array([500,550,600,650.0]) + 273.15

  for T in Ts:
    print(T)
    print("Aging")
    drv = drivers.Driver_sd(model, verbose = True, T_init = T)
    for t in np.logspace(-1,np.log10(aging),nage):
      drv.strain_step(np.zeros((6,)), t, T) 
    
    print("Pulling")
    res = drivers.uniaxial_test(model, 8.33e-5, T = T, verbose = True,
        history = drv.stored_int[-1])
  
    plt.plot(res['strain'], res['stress'], label = "%3.0f C" % (T-273.15))

  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  plt.savefig("tension-aged.png")

