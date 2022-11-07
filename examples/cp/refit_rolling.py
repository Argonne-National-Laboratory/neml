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
from neml import elasticity, parse

import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Number of crystals and number of threads
    N = 100
    nthreads = 1

    # Temperature
    T = 298.0

    # Sets up the lattice crystallography
    a = 2.9511 * 0.1  # nm
    c = 4.68433 * 0.1  # nm
    lattice = crystallography.HCPLattice(a, c)

    # Strain direction, rate, and number of steps
    L = np.array([[0.0, 0, 0], [0, 1.0, 0], [0, 0, -1.0]])
    erate = 1.0e-3
    steps = 100
    emax = 2.0

    # Sets up the actual deformation tensor
    L *= erate
    dt = emax / steps / erate
    # Randomly selects initial orientations
    orientations = rotations.random_orientations(N)

    model = parse.parse_xml("refit.xml", "ti_cp")
    # Sets up the poly crystal model
    pmodel = polycrystal.TaylorModel(model, orientations, nthreads=nthreads)

    # Runs the exmaple
    h_n = pmodel.init_store()

    d_inc = nemlmath.sym(0.5 * (L + L.T))
    w_inc = nemlmath.skew(0.5 * (L - L.T))

    s_n = np.zeros((6,))
    d_n = np.zeros((6,))
    w_n = np.zeros((3,))

    u_n = 0.0
    p_n = 0.0

    for i in range(steps):
        print(i)
        d_np1 = d_n + d_inc * dt
        w_np1 = w_n + w_inc * dt
        s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = pmodel.update_ld_inc(
            d_np1, d_n, w_np1, w_n, T, T, dt, 0, s_n, h_n, u_n, p_n
        )

        d_n = np.copy(d_np1)
        w_n = np.copy(w_np1)

        s_n = np.copy(s_np1)
        h_n = np.copy(h_np1)

        u_n = u_np1
        p_n = p_np1

    # Plots a second, as-rolled basal pole figure
    polefigures.pole_figure_discrete(pmodel.orientations(h_np1), [0, 0, 0, 1], lattice)
    plt.title("Final, <0001>")
    plt.show()
    plt.close()

    polefigures.pole_figure_discrete(pmodel.orientations(h_np1), [1, 0, -1, 0], lattice)
    plt.title("Final, <1010>")
    plt.show()
    plt.close()

    polefigures.pole_figure_discrete(pmodel.orientations(h_np1), [1, 0, -1, 1], lattice)
    plt.title("Final, <1011>")
    plt.show()
    plt.close()
