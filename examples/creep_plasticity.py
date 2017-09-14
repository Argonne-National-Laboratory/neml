#!/usr/bin/env python

import sys
sys.path.append('..')
sys.path.append('../test')

from neml import solvers, neml, elasticity, drivers, surfaces, hardening, ri_flow, creep

from common import differentiate

import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la

if __name__ == "__main__":
    A = 1.85e-10
    n = 2.5

    smodel = creep.PowerLawCreep(A, n)
    cmodel = creep.J2CreepModel(smodel)

    E = 150000.0
    nu = 0.3
    sY = 200.0

    youngs = elasticity.YoungsModulus(E)
    poisson = elasticity.PoissonsRatio(nu)
    elastic = elasticity.IsotropicLinearElasticModel(youngs, 
        poisson)
    surface = surfaces.IsoJ2()

    pmodel = neml.SmallStrainPerfectPlasticity(elastic, surface, sY)

    model = neml.SmallStrainCreepPlasticity(pmodel, cmodel)
    
    """
    res = drivers.creep(model, 200.0, 1.0e2, 1000.0, verbose = True)
    plt.plot(res['time'], res['strain'])
    plt.show()

    res = drivers.strain_cyclic(model, 0.01, -1, 1.0e-4, 15, verbose = True,
        nsteps = 25, hold_time = [10.0,0])
    plt.plot(res['strain'], res['stress'])
    plt.show()
    
    res = drivers.uniaxial_test(model, 1.0e-6, sdir = np.array([1,-1,0.5,-0.5,0.1,-0.25]), 
        verbose = True)
    plt.plot(res['strain'], res['stress'])
    plt.show()
    """
    
    tfinal = 100
    T = 300.0
    nsteps = 10

    t_n = 0.0
    strain_n = np.zeros((6,))
    stress_n = np.zeros((6,))
    hist_n = model.init_store()

    efinal = np.array([0,0,0,0.0,0.1,0.0])

    u_n = 0.0
    p_n = 0.0

    for i,m in enumerate(np.linspace(0,1,nsteps)):
      t_np1 = tfinal * m
      strain_np1 = efinal * m
      
      stress_np1, hist_np1, A_np1, u_np1, p_np1 = model.update_sd(
          strain_np1, strain_n, T, T, t_np1, t_n, stress_n, hist_n,
          u_n, p_n)
      dfn = lambda e: model.update_sd(e,
          strain_n, T, T, t_np1, t_n, stress_n, hist_n, u_n, p_n)[0]

      print("HMM")
      print(A_np1)
      #print(num_A)

      num_A = differentiate(dfn, strain_np1, eps = 1.0e-6)
      
      strain_n = strain_np1
      stress_n = stress_np1
      hist_n = hist_np1
      u_n = u_np1
      p_n = p_np1
