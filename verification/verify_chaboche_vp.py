#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import solvers, neml, elasticity, drivers, surfaces, hardening, visco_flow, general_flow

import matplotlib.pyplot as plt
import numpy as np

def koo():
  # Data from Koo & Kwon (2011)
  E = 157.0e3
  nu = 0.27

  sY = 0.0
  b = 3.0
  Q = -80.0 

  C1 = 135.0e3
  C2 = 123.0e3
  C3 = 4.0e3

  g1 = 100e3
  g2 = 0.85e3
  g3 = 1.0

  eta = 701.0

  n = 10.5

  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  shear = elasticity.ShearModulus(mu)
  bulk = elasticity.BulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

  surface = surfaces.IsoKinJ2()
  iso = hardening.VoceIsotropicHardeningRule(sY, Q, b)
  cs = np.array([C1, C2, C3])
  gs = np.array([g1, g2, g3])
  hmodel = hardening.Chaboche(iso, cs, gs)

  fluidity = visco_flow.ConstantFluidity(eta)

  vmodel = visco_flow.ChabocheFlowRule(
      surface, hmodel, fluidity, n)

  flow = general_flow.TVPFlowRule(elastic, vmodel)

  model = neml.GeneralIntegrator(flow)


  e500 = np.array([0.0058793164, 0.0052876081, 0.004710156, 0.0042916722, 0.0038731605, 0.0032237573, 0.0027900399, 0.0024864991, 0.0017489338, 0.000635346, -0.0006668255, -0.0017811254, -0.0034891722, -0.0044012326, -0.0059358454, -0.005661599, -0.0051996206, -0.0048027165, -0.0044347716, -0.0039945684, -0.0034388636, -0.0032440594, -0.0026585856, -0.002036857, -0.001161945, -0.0003736313, 0.0008344368, 0.0022454705, 0.0035554683, 0.0052854997])
  s500 = np.array([349.6965238343, 255.0932943767, 171.2928564634, 102.3697757811, 34.7967386021, -70.6166129529, -103.0909625627, -128.7931029969, -172.1191145041, -237.1093535238, -285.2558144812, -315.8199441656, -342.4341275638, -350.0133084444, -361.0726693423, -319.1748939033, -252.9445696747, -190.7753717856, -128.6110609317, -59.6843149729, 25.4625011672, 57.8965327357, 103.8969743143, 147.2034376803, 200.6779888803, 240.6374439409, 288.0930002819, 323.4323743399, 343.229143486, 357.0215788789]) 

  res = drivers.strain_cyclic(model, 0.006, -1.0, 1.0e-4, 130,
      verbose = False, nsteps = 50)
  plt.plot(res['strain'][-100:], res['stress'][-100:], 'k-')
  plt.plot(e500, s500, 'kx')
  plt.xlabel("Strain (-/-)")
  plt.ylabel("Stress (MPa)")
  plt.show()

  n500 = np.array([1.8000456006, 7.2646606231, 17.3775914892, 25.6289731945, 35.040864809, 43.9864865417, 54.6731673353, 69.3076477381, 83.9405590891, 98.1077889527, 112.0424068927, 126.7891398736])
  ns500 = sn500 = np.array([437.8476941068, 428.2481618804, 413.6055992266, 403.8518586317, 394.8643791744, 388.2393032756, 382.3835530687, 376.4641827161, 372.5718311079, 370.2817510291, 368.4971998147, 367.5615301197])

  plt.plot(res['cycles'], res['max'], 'k-')
  plt.plot(n500, ns500, 'kx')
  plt.xlabel("Cycle")
  plt.ylabel("Stress (MPa)")
  plt.show()


def uniaxial():
  # Data from "Numerical modeling of elasto-viscoplastic Chaboche constitutive..."
  # by A. Ambroziak
  E = 159000.0
  nu = 0.3

  k = 514.21
  b = 60.0
  R1 = -194.39
  a = 170000.0
  c = 500.0
  n = 4.0
  eta = 1023.5


  # Translate
  Q = R1
  b = b
  
  C1 = a
  g1 = c

  eta = eta
  sY = k

 
  mu = E / (2 * (1.0 + nu))
  K = E / (3 * (1 - 2 * nu))

  shear = elasticity.ShearModulus(mu)
  bulk = elasticity.BulkModulus(K)
  elastic = elasticity.IsotropicLinearElasticModel(shear, bulk)

  surface = surfaces.IsoKinJ2()
  iso = hardening.VoceIsotropicHardeningRule(sY, Q, b)
  cs = np.array([C1])
  gs = np.array([g1])
  hmodel = hardening.Chaboche(iso, cs, gs)

  fluidity = visco_flow.ConstantFluidity(eta)

  vmodel = visco_flow.ChabocheFlowRule(
      surface, hmodel, fluidity, n)

  flow = general_flow.TVPFlowRule(elastic, vmodel)

  model = neml.GeneralIntegrator(flow)

  erates = [1.0e-7, 1.0e-2, 1.0e-1]
 
  e_7 = np.array([1.03E-007, 0.000977076, 0.0019927018, 0.0029695719, 0.0038278073, 0.0048024128, 0.0058155682, 0.0068186871, 0.0078116926, 0.0088238185, 0.0098260368, 0.0108280492, 0.0118396862, 0.0128318167, 0.0138433765, 0.0148451316, 0.0158469381, 0.0168487447, 0.017840798, 0.0188523063, 0.019844411])
  s_7 = np.array([1.9738128815, 164.2845378653, 316.6592895303, 475.0018013938, 559.9786920852, 631.0205153016, 688.1581520788, 733.3950096246, 763.7549280987, 781.2104336727, 791.7254264156, 794.3039929179, 792.9104861704, 789.5405931218, 785.170926534, 777.8289602355, 772.4711004972, 767.1132407589, 760.7671878699, 754.413414722, 750.0514683933])
 
  e_2 = np.array([-1.94E-005, 0.000976973, 0.0019830515, 0.0029890784, 0.0037900553, 0.0048443338, 0.0058197114, 0.0068040701, 0.0078172512, 0.0088300977, 0.0098231803, 0.0108257074, 0.0118377047, 0.0128397943, 0.0138612103, 0.0148630941, 0.0158551988, 0.0168667329, 0.017878267, 0.0188703202, 0.0198623735])
  s_2 = np.array([-0.0025734197, 160.316324745, 319.6393095, 476.9781876949, 604.6365302068, 747.075308553, 847.8787301718, 919.9087465388, 978.0384365961, 1023.2714340124, 1056.6075123267, 1079.0271444306, 1091.5223836042, 1097.0771099468, 1097.66384963, 1095.2821497319, 1090.9202034031, 1085.5584835353, 1080.1967636674, 1073.8507107785, 1067.5046578896])

  e_1 = np.array([-1.93E-005, 0.0009672712, 0.0019928047, 0.0029793251, 0.0037900295, 0.0047863804, 0.0058017489, 0.0068163453, 0.0078302727, 0.0088243075, 0.0098374371, 0.0108404531, 0.0118432376, 0.0128554407, 0.0138577362, 0.014850124, 0.01587154, 0.0168734495, 0.0178753075, 0.0188770883, 0.0198789206])
  s_1 = np.array([2.9735864206, 161.3122381546, 320.6275026506, 475.9899945444, 603.6444769267, 762.9713218113, 905.4255406755, 1018.1181611373, 1105.0173963169, 1175.0594459942, 1231.2050294914, 1272.4736739169, 1304.8138388216, 1325.2455042359, 1338.7366568191, 1345.2872965712, 1345.8740362543, 1344.4843896363, 1341.1106364582, 1334.7607234397, 1330.3949169815])

  plt.plot(e_7, s_7, 'kx')
  plt.plot(e_2, s_2, 'rx')
  plt.plot(e_1, s_1, 'bx')

  res = drivers.uniaxial_test(model, erates[0], emax = 0.02)
  plt.plot(res['strain'], res['stress'], 'k-')

  res = drivers.uniaxial_test(model, erates[1], emax = 0.02)
  plt.plot(res['strain'], res['stress'], 'r-')

  res = drivers.uniaxial_test(model, erates[2], emax = 0.02)
  plt.plot(res['strain'], res['stress'], 'b-')

  plt.xlim([0,0.02])
  plt.ylim([0,1500])

  plt.show()

if __name__ == "__main__":
  uniaxial()
  koo()
