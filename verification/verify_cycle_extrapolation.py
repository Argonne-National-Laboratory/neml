#!/usr/bin/env python3

import sys
sys.path.append('..')

import numpy as np
from neml import interpolate,solvers, drivers,models, elasticity, ri_flow, hardening, surfaces, visco_flow, general_flow, creep, damage

import matplotlib.pyplot as plt
import time

if __name__ == "__main__":

    E = 150000
    nu= 0.3

    n = 9.282
    eta = 228.09
    sY = 7.485

    Q = 153.9
    b = 1.9878

    C1 = 15118
    g1 = 309.63

    C2 = 2782.8
    g2 = 18.239

    surface = surfaces.IsoKinJ2()
    iso = hardening.VoceIsotropicHardeningRule(sY, Q, b)
    cs = [C1, C2]
    gs = np.array([g1, g2])
    gmodels = [hardening.ConstantGamma(g) for g in gs]
    hmodel = hardening.Chaboche(iso, cs, gmodels, [0.0] * len(cs),[1.0] * len(cs))

    fluidity = visco_flow.ConstantFluidity(eta)

    vmodel = visco_flow.ChabocheFlowRule(surface, hmodel, fluidity, n)

    mu = E/(2*(1+nu))
    K = E/(3*(1-2*nu))

    elastic = elasticity.IsotropicLinearElasticModel(mu, "shear", K, "bulk")
    fmodel = general_flow.TVPFlowRule(elastic, vmodel)
    bmodel = models.GeneralIntegrator(elastic, fmodel, verbose = False)
    A          = 1e4
    xi         = 4.0
    phi        = 2.0

    model = damage.ModularCreepDamageModel_sd(elastic, A, xi, phi,damage.VonMisesEffectiveStress(),bmodel)

    emax    = np.array([ 2.0  , 2.0  , 1.0   , 2.0  ])/200.0
    R       = np.array([ -1.0 , -1.0 , -1.0  , -1.0 ])
    erate   = np.array([ 0.004, 0.004, 0.004 , 0.004])
    hold_t  = np.array([ 0.0  , 6.0  , 30.0  , 6.0  ])*60.0
    hold_c  = np.array([ 0.0  , 0.0  ,  0.0  , 6.0  ])*60.0
    ncycles = np.array([ 500  , 500  , 200   , 120  ])

    # time_act = []
    # time_ext = []
    # for r in range(len(emax)):
    #     start = time.time()
    #     res = drivers.strain_cyclic_extrapolated(model, emax[r], R[r], erate[r], int(ncycles[r]),hold_time=[hold_t[r],hold_c[r]],
    #                          min_cycle=3, unit_extrapolate = 10,allowable_jump_stress=10.0,jump_del_N=10,check_dmg=True)
    #     end = time.time()
    #     time_ext.append(end - start)
    #     plt.plot(res['cycles'],res['max'],'+',label='ext run {}'.format(r+1))
    #     plt.plot(res['cycles'],res['min'],'+')
    #     start = time.time()
    #     res = drivers.strain_cyclic(model, emax[r], R[r], erate[r], int(ncycles[r]),hold_time=[hold_t[r],hold_c[r]],check_dmg=True)
    #     end = time.time()
    #     time_act.append(end - start)
    #     plt.plot(res['cycles'],res['max'],label='act run {}'.format(r+1))
    #     plt.plot(res['cycles'],res['min'])
    # plt.legend(loc='best')
    # plt.xlabel('cycles')
    # plt.ylabel('stress')
    # plt.title('Extrapolation comparison with damage')
    # plt.show()
    # plt.scatter(np.arange(len(emax)),time_act,label='act')
    # plt.scatter(np.arange(len(emax)),time_ext,marker='^',label='ext')
    # plt.yscale('log')
    # plt.legend(loc='best')
    # plt.ylabel('time')
    # plt.xlabel('runs')
    # plt.title('Time comparison')
    # plt.show()

    ## no damage
    time_act = []
    time_ext = []
    for r in range(len(emax)):
        start = time.time()
        res = drivers.strain_cyclic_extrapolated(bmodel, emax[r], R[r], erate[r], int(ncycles[r]),hold_time=[hold_t[r],hold_c[r]],
                                                   min_cycle=3, unit_extrapolate = 10,allowable_jump_stress=10.0,jump_del_N=10)
        end = time.time()
        time_ext.append(end - start)
        plt.plot(res['cycles'],res['max'],'+',label='ext run {}'.format(r+1))
        plt.plot(res['cycles'],res['min'],'+')
        start = time.time()
        res = drivers.strain_cyclic(bmodel, emax[r], R[r], erate[r], int(ncycles[r]),hold_time=[hold_t[r],hold_c[r]])
        end = time.time()
        time_act.append(end - start)
        plt.plot(res['cycles'],res['max'],label='act run {}'.format(r+1))
        plt.plot(res['cycles'],res['min'])
    plt.legend(loc='best')
    plt.xlabel('cycles')
    plt.ylabel('stress')
    plt.title('Extrapolation comparison no damage')
    plt.show()
    plt.scatter(np.arange(len(emax)),time_act,label='act')
    plt.scatter(np.arange(len(emax)),time_ext,marker='^',label='ext')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.ylabel('time')
    plt.xlabel('runs')
    plt.title('Time comparison')
    plt.show()
