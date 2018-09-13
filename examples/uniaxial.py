#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow, uniaxial

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
  E = 200000.0
  nu = 0.27

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  s0 = 300.0
  Kp = 0.0
  c = [30000.0]
  r = [60.0]
  
  elastic = elasticity.IsotropicLinearElasticModel(mu, "shear", K, "bulk")
  surface = surfaces.IsoKinJ2()
  iso = hardening.LinearIsotropicHardeningRule(s0, Kp)
  gmodels = [hardening.ConstantGamma(g) for g in r]
  hrule = hardening.Chaboche(iso, c, gmodels, [0.0] * len(c), 
      [1.0] * len(c))

  flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hrule)
  model = models.SmallStrainRateIndependentPlasticity(elastic, flow, verbose = False,
      check_kt = False)
  umodel = uniaxial.UniaxialModel(model, verbose = True)

  hn = umodel.init_store()
  en = 0.0
  sn = 0.0
  un = 0.0
  pn = 0.0
  tn = 0.0
  Tn = 0.0

  dt = 1.0
  dT = 0.0

  es = [en]
  ss = [sn]

  n = 100
  emax = 0.1

  for enp1 in np.linspace(0, emax, n+1)[1:]:
    Tnp1 = Tn + dT
    tnp1 = tn + dt
    snp1, hnp1, Anp1, unp1, pnp1 = umodel.update(enp1, en, Tnp1, Tn, 
        tnp1, tn, sn, hn, un, pn)
    es.append(enp1)
    ss.append(snp1)

    sn = snp1
    hn = np.copy(hnp1)
    en = enp1
    Tn = Tnp1
    tn = tnp1
    un = unp1
    pn = pnp1

  plt.plot(es, ss, 'k-')

  res = drivers.uniaxial_test(model, emax / n / dt, T = Tn, emax = emax)
  
  plt.plot(res['strain'], res['stress'], 'r-')
  plt.show()
    
