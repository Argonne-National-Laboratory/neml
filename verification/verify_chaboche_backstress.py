#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, models, elasticity, drivers, surfaces, hardening, ri_flow

import matplotlib.pyplot as plt
import numpy as np

def verify_Q():
  E = 30000.0
  nu = 0.3

  sy = 100.0
  Q = 50.0
  b = 100.0
  
  C = 0.0
  g = 0.0

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  shear = elasticity.ShearModulus(mu)
  bulk = elasticity.BulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)
  surface = surfaces.IsoKinJ2()
  iso = hardening.VoceIsotropicHardeningRule(sy, Q, b)
  gmodels = [hardening.ConstantGamma(g)]
  hrule = hardening.Chaboche(iso, [C], gmodels)

  flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hrule)
  model = models.SmallStrainRateIndependentPlasticity(elastic, flow, verbose = False,
      check_kt = False)

  res = drivers.uniaxial_test(model, 1.0e-2, emax = 0.2)
  stress = res['stress']

  print("Q: %f / %f" % (Q, stress[-1] - sy))


def verify_Cg():
  E = 30000.0
  nu = 0.3

  sy = 100.0
  Q = 0.0
  b = 0.0
  
  C = 1000.0
  g = 10.0

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  shear = elasticity.ShearModulus(mu)
  bulk = elasticity.BulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)
  surface = surfaces.IsoKinJ2()
  iso = hardening.VoceIsotropicHardeningRule(sy, Q, b)
  hrule = hardening.Chaboche(iso, [C], [hardening.ConstantGamma(g)])

  flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hrule)
  model = models.SmallStrainRateIndependentPlasticity(elastic, flow, verbose = False,
      check_kt = False)

  res = drivers.strain_cyclic(model, 0.4, -1.0, 1.0e-4, 1)
  strain = res['strain']
  stress = res['stress']
  
  mv = np.max(np.abs(stress))
  hu = mv - sy

  print("C/y: %f / %f" % ((C/g), hu))
   
def verify_warp3d():
  E = 30000.0
  nu = 0.3

  sy = 100.0
  Q = 50.0
  b = 100.0
  
  C = 1000.0
  g = 10.0

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  shear = elasticity.ShearModulus(mu)
  bulk = elasticity.BulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)
  surface = surfaces.IsoKinJ2()
  iso = hardening.VoceIsotropicHardeningRule(sy, Q, b)
  hrule = hardening.Chaboche(iso, [C], [hardening.ConstantGamma(g)])

  flow = ri_flow.RateIndependentNonAssociativeHardening(surface, hrule)
  model = models.SmallStrainRateIndependentPlasticity(elastic, flow, verbose = False,
      check_kt = False)

  res = drivers.strain_cyclic(model, 0.0075, -1.0, 1.0e-4, 50)
  strain = res['strain']
  stress = res['stress']

  data_warp = np.load('data_fa_warp.npy')

  plt.plot(strain, stress, 'k-')
  plt.plot(data_warp[0], data_warp[1], 'r-')
  plt.show()

if __name__ == "__main__":
  verify_Q()
  
  verify_Cg()

  verify_warp3d()
