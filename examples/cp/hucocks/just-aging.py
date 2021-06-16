#!/usr/bin/env python3

from generate_model import *

from neml import drivers, history

import matplotlib.pyplot as plt

fs = 0.1
rs = 1.0e-9
Ns = 1.0e12

L = crystallography.CubicLattice(1.0)
L.add_slip_system([1,1,0],[1,1,1])

def integrate_case(model, times, T):
  driver = drivers.Driver_sd(model, verbose = True, T_init = T)
  
  for t in times[1:]:
      driver.strain_step(np.zeros((6,)), t, T)
  
  store = np.array(driver.stored_int)

  return store, store[:,-6:-3], store[:,-3:]

if __name__ == "__main__":
  model, hmodel = make_singlecrystal(verbose = True, return_hardening = True)

  times = np.insert(np.logspace(0,7,100) * 3600.0, 0, [0])
  times = np.insert(times, 1, [1.0e2])
 
  Ts = np.array([500,550,600,650.0]) + 273.15

  for T in Ts:
    print(T)
    h, carbide, laves = integrate_case(model, times, T)
    
    hours = times / 3600.0
    f_laves = laves[:,0] * fs
    r_laves = laves[:,1] * rs
    N_laves = laves[:,2] * Ns
    A_laves = 2.0 * r_laves * N_laves
    
    f_carbide = carbide[:,0] * fs
    r_carbide = carbide[:,1] * rs
    N_carbide = carbide[:,2] * Ns
    A_carbide = 2.0 * r_carbide * N_carbide

    taus = []
    for hi in h:
      blank = history.History()
      hist = history.History()
      model.populate_history(hist)
      hist.set_data(hi)
      taus.append(hmodel.hist_to_tau(0, 0, hist, L, T, blank))
  
    plt.figure()
    plt.semilogx(hours, f_laves, label = "Laves")
    plt.semilogx(hours, f_carbide, label = "Carbide")
    plt.legend(loc='best')
    plt.xlabel("Time (hours)")
    plt.ylabel("Volume fraction")
    plt.xlim([1e0,1e7])
    plt.show()

    plt.figure()
    plt.semilogx(hours, r_laves * 1.0e9, label = "Laves")
    plt.semilogx(hours, r_carbide * 1.0e9, label = "Carbide")
    plt.legend(loc='best')
    plt.xlabel("Time (hours)")
    plt.ylabel("Radius (nm)")
    plt.xlim([1e0,1e7])
    plt.show()

    plt.figure()
    plt.semilogx(hours, A_laves, label = "Laves")
    plt.semilogx(hours, A_carbide, label = "Carbide")
    plt.legend(loc='best')
    plt.xlabel("Time (hours)")
    plt.ylabel("Area density (1/m^2)")
    plt.xlim([1e0,1e7])
    plt.show()

    plt.figure()
    plt.semilogx(hours, taus)
    plt.xlabel("Time (hours)")
    plt.ylabel("Slip resistance (MPa)")
    plt.xlim([1e0,1e7])
    plt.show()
