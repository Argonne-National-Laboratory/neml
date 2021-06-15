#!/usr/bin/env python3

from generate_model import *

from neml import drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 1
  nthreads = 2
  model = make_model(N, nthreads = nthreads)
  
  res = drivers.uniaxial_test(model, 8.33e-5, T = 600+273.15)
  
  plt.plot(res['strain'], res['stress'])
  plt.show()

