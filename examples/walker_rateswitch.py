#!/usr/bin/env python3

import sys
sys.path.append('..')

from neml import models, elasticity, walker, general_flow, drivers, surfaces, ri_flow, hardening

import matplotlib.pyplot as plt

if __name__ == "__main__":
  E = 140000.0
  nu = 0.33

  eps0 = 1.0e-4
  eps_ref = 1.0e-10
  D = 100.0
  n = 5.62
  s0 = 150.0
  K = E / 50.0

  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  vmodel = walker.TestFlowRule(eps0, D, n, s0, K)
  frule = general_flow.TVPFlowRule(emodel, vmodel)
  model1 = models.GeneralIntegrator(emodel, frule)

  frule2 = walker.WalkerKremplSwitchRule(emodel, vmodel, 0.0, eps_ref)
  model2 = models.GeneralIntegrator(emodel, frule2)

  frule3 = walker.WalkerKremplSwitchRule(emodel, vmodel, 0.99, eps_ref)
  model3 = models.GeneralIntegrator(emodel, frule3)

  rates = [1.0e-6,1.0e-4,1.0e-2,1.0]

  surface = surfaces.IsoJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, K)
  flow = ri_flow.RateIndependentAssociativeFlow(surface, iso)
  model4 = models.SmallStrainRateIndependentPlasticity(emodel, flow)

  for r in rates:
    res1 = drivers.uniaxial_test(model1, r)
    l, = plt.plot(res1['strain'], res1['stress'])
    res2 = drivers.uniaxial_test(model2, r)
    plt.plot(res2['strain'], res2['stress'], color = l.get_color(), ls = '--', lw = 2)
    res3 = drivers.uniaxial_test(model3, r)
    plt.plot(res3['strain'], res3['stress'], color = l.get_color(), ls = ':', lw = 3)
    res4 = drivers.uniaxial_test(model4, r)
    plt.plot(res4['strain'], res4['stress'], color = 'k', ls = '-')
    """
    print("Verification")
    print("E = %f versus %f" % ((res['stress'][1] - res['stress'][0]) / (res['strain'][1] - res['strain'][0]), E))
    print("K = %f versus %f" % ((res['stress'][-2] - res['stress'][-1]) / (res['strain'][-2] - res['strain'][-1]), K))
    print("eps = %f versus %f" % (res['history'][0], res['strain'][-1] - res['stress'][-1] / E))
    """

  plt.show()
