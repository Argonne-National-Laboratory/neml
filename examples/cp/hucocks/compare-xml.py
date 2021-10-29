#!/usr/bin/env python3

from generate_model import *

from neml import drivers, models, parse

import matplotlib.pyplot as plt

if __name__ == "__main__":
  model1 = make_singlecrystal()
  model2 = parse.parse_xml("model.xml", "hucocks")

  Ts = np.array([500,550,600,650.0]) + 273.15
  ss = np.array([25.0,25.0,25.0,25.0])
  time = 50*3600.0
  srate = 1.0

  for T, s in zip(Ts,ss):
    print(T)
    res1 = drivers.creep(model1, s, srate, time, T = T, 
        nsteps = 50, nsteps_up = 50, logspace = False, verbose = True)
    res2 = drivers.creep(model2, s, srate, time, T = T, 
        nsteps = 50, nsteps_up = 50, logspace = False, verbose = True)
    
    plt.figure()
    plt.plot(res1['time']/3600.0, res1['strain'])
    plt.plot(res2['time']/3600.0, res2['strain'])
    plt.show()

  
