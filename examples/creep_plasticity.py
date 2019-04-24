#!/usr/bin/env python3

import sys
sys.path.append('..')
sys.path.append('../test')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow, creep

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

    elastic = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, 
        "poissons")
    surface = surfaces.IsoJ2()

    pmodel = models.SmallStrainPerfectPlasticity(elastic, surface, sY)

    model = models.SmallStrainCreepPlasticity(elastic, pmodel, cmodel)
    
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
