#!/usr/bin/env python3

import sys

sys.path.append("..")

from neml import (
    solvers,
    models,
    elasticity,
    drivers,
    surfaces,
    hardening,
    ri_flow,
    interpolate,
)

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    E = 200000.0
    nu = 0.27

    temperatures = [0.0, 100.0, 200.0, 300.0]
    points = [0.0, 0.01, 0.02, 0.045]

    values = [
        [50.0, 75.0, 100.0, 120.0],
        [100.0, 150.0, 160.0, 170.0],
        [200.0, 250.0, 260.0, 270.0],
        [300.0, 350.0, 360.0, 370.0],
    ]

    flow_curves = [interpolate.PiecewiseLinearInterpolate(points, v) for v in values]

    elastic = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
    surface = surfaces.IsoJ2()
    iso = hardening.TemperatureDependentInterpolatedIsotropicHardeningRule(
        temperatures, flow_curves
    )
    flow = ri_flow.RateIndependentAssociativeFlow(surface, iso)
    model = models.SmallStrainRateIndependentPlasticity(elastic, flow)

    # Play around with test_temperature to make sure you like the behavior, try values above and below the interpolation range
    test_temperature = 150.0
    res = drivers.uniaxial_test(model, 1.0e-4, emax=0.045, T=test_temperature)

    plt.plot(res["strain"], res["stress"], "k-")

    plot_strains = np.linspace(0.0, 0.045, 100)
    for i, T in enumerate(temperatures):
        print(i, T)
        plt.plot(
            plot_strains + values[i][0] / E,
            [flow_curves[i](e) for e in plot_strains],
            "r-",
        )

    plt.xlim([0, 0.05])
    plt.show()
