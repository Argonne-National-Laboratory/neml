#!/usr/bin/env python

import numpy as np
import scipy.interpolate as inter
import scipy.optimize as opt
import os.path

import sys
sys.path.append('../..')

from neml import solvers, neml, interpolate, elasticity, drivers, surfaces, hardening, ri_flow, visco_flow, general_flow

def make_yaguchi_model():
  vmodel = visco_flow.YaguchiGr91FlowRule()

  E_poly = np.array([-3.077e-1, 3.003e2, 1.269e5])
  nu = 0.3

  mu_poly = list(E_poly/(2*(1+nu)))
  K_poly = list(E_poly/(3*(1-2*nu)))

  shear = elasticity.ShearModulus(interpolate.PolynomialInterpolate(mu_poly))
  bulk = elasticity.BulkModulus(interpolate.PolynomialInterpolate(K_poly))
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

  flow = general_flow.TVPFlowRule(elastic, vmodel)

  return neml.GeneralIntegrator(flow)

def make_koo_model(temp):
  """
    Return a fresh model at the given temperature
  """
  if temp == "500":
    E = 157.0e3

    C1 = 135.0e3
    C2 = 123.2e3
    C3 = 4.0e3
    g1 = 1.0e5
    g2 = 850.0
    g3 = 1.0

    b = 3.0
    Q = -80.0
    K = 701.0

    n = 10.5

  elif temp == "550":
    E = 155.0e3

    C1 = 135.0e3
    C2 = 107.2e3
    C3 = 6.0e3
    g1 = 1.0e5
    g2 = 1000.0
    g3 = 1.0

    b = 4.0
    Q = -100.0
    K = 1265.0

    n = 6.0

  elif temp == "600":
    E = 150.0e3

    C1 = 135.0e3
    C2 = 61.0e3
    C3 = 11.0e3
    g1 = 5.0e4
    g2 = 1100.0
    g3 = 1.0

    b = 5.0
    Q = -110.0
    K = 3075.0

    n = 3.6

  else:
    raise ValueError("Unknown temperature %s." % temp)

  surface = surfaces.IsoKinJ2()
  iso = hardening.VoceIsotropicHardeningRule(0.0, Q, b)
  cs = [C1, C2, C3]
  gs = [g1, g2, g3]
  gmodels = [hardening.ConstantGamma(g) for g in gs]
  hmodel = hardening.Chaboche(iso, cs, gmodels)

  fluidity = visco_flow.ConstantFluidity(K)

  vmodel = visco_flow.ChabocheFlowRule(surface, hmodel, fluidity, n)

  nu = 0.3

  mu = E/(2*(1+nu))
  K = E/(3*(1-2*nu))

  shear = elasticity.ShearModulus(mu)
  bulk = elasticity.BulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

  flow = general_flow.TVPFlowRule(elastic, vmodel)

  return neml.GeneralIntegrator(flow, verbose = False)
