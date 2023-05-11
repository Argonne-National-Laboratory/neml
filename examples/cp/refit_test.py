#!/usr/bin/env python3

import sys

sys.path.append("../..")

import numpy as np

from neml.cp import (
    polycrystal,
    crystallography,
    slipharden,
    sliprules,
    inelasticity,
    kinematics,
    singlecrystal,
    polefigures,
    postprocessors,
)
from neml.math import rotations, tensors, nemlmath, matrix
from neml import elasticity, parse, drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Number of crystals and number of threads
    N = 1
    nthreads = 1

    # Randomly selects initial orientations
    orientations = rotations.random_orientations(N)

    model = parse.parse_xml("refit.xml", "ti_cp")
    # Sets up the poly crystal model
    pmodel = polycrystal.TaylorModel(model, orientations, nthreads=nthreads)

    res = drivers.uniaxial_test(
        model,
        erate=1.0e-2,
        emax=0.2,
        sdir=np.array([0, 0, 1, 0, 0, 0]),
        T=298.15,
        verbose=True,
        full_results=False,
    )

    plt.plot(res["strain"], res["stress"], lw=2, alpha=0.5)
    plt.show()