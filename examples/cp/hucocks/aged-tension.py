#!/usr/bin/env python3

from generate_model import *

from neml import drivers

import matplotlib.pyplot as plt

def integrate_case(model, time, T):
  driver = drivers.Driver_sd(model, verbose = True, T_init = T)

  times = np.insert(np.logspace(0,np.log10(time),100) * 3600.0, 0, [0])
  times = np.insert(times, 1, [1.0e2])

  for t in times[1:]:
      driver.strain_step(np.zeros((6,)), t, T)
  
  store = np.array(driver.stored_int)

  return store

if __name__ == "__main__":
  N = 100
  nthreads = 3
  model = make_model(N, nthreads = nthreads)
 
  Ts = np.array([500,550,600,650.0]) + 273.15
  
  for T in Ts:
    print(T)
    hist = integrate_case(model, 100000, T)
    res = drivers.uniaxial_test(model, 8.33e-5, T = T, verbose = True, 
        nsteps = 100, emax = 0.030, history = hist)
    np.savetxt("aged-%i.txt" % T, np.stack((res['strain'], res['stress'])).T)
  


