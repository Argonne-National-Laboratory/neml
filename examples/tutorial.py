#!/usr/bin/env python3

from neml import solvers, models, elasticity, drivers, surfaces, hardening, visco_flow, general_flow, parse

import matplotlib.pyplot as plt

E = 150000.0
nu = 0.3
elastic = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, 
    "poissons")

surface = surfaces.IsoJ2()

sY = 100.0
H = 2500.0
iso = hardening.LinearIsotropicHardeningRule(sY, H)

eta = 100.0
n = 5.0
gpower = visco_flow.GPowerLaw(n, eta)
vflow = visco_flow.PerzynaFlowRule(surface, iso, gpower)

integrator = general_flow.TVPFlowRule(elastic, vflow)

model = models.GeneralIntegrator(elastic, integrator)

erate = 1.0e-4
res = drivers.uniaxial_test(model, erate)

plt.figure()
plt.plot(res['strain'], res['stress'], 'k-')
plt.xlabel("Strain (mm/mm)")
plt.ylabel("Stress (MPa)")
plt.show()

model2 = parse.parse_xml("tutorial.xml", "tutorial_model")

res2 = drivers.uniaxial_test(model2, erate)

plt.figure()
plt.plot(res2['strain'], res2['stress'], 'k-')
plt.xlabel("Strain (mm/mm)")
plt.ylabel("Stress (MPa)")
plt.show()
