#!/usr/bin/env python3

from generate_model_Ti import *

from neml import drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 1000
  nthreads = 1
  strain_rate = 1.0e-4
  model = make_model(strain_rate, N, nthreads = nthreads)
 
  T = 298.0
  res = drivers.uniaxial_test(model, strain_rate, T = T, verbose = True)
  
  plt.plot(res['strain'], res['stress'], label = "%3.0f C" % (T-273.15))

  """
  Ts = np.array([500,550,600,650.0]) + 273.15
  for T in Ts:
    res = drivers.uniaxial_test(model, strain_rate, T = T, verbose = True)
  
    plt.plot(res['strain'], res['stress'], label = "%3.0f C" % (T-273.15))
  """
  
  plt.legend(loc='best')
  plt.xlabel("Strain (mm/mm)")
  plt.ylabel("Stress (MPa)")
  # plt.savefig("tension-Ti.png")
  plt.show()
  plt.close()

