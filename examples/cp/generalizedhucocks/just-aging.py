#!/usr/bin/env python

import sys

sys.path.append("../../..")

from neml import parse
from neml import drivers, history
from pathlib import Path

from generate_model import *
import matplotlib.pyplot as plt

import pdb

rs = 1.0e-9
Ns = 1.0e12

L = crystallography.CubicLattice(1.0)
L.add_slip_system([1, 1, 0], [1, 1, 1])

if __name__ == "__main__":
    # smodel = parse.parse_xml(str(Path(__file__).parent / "model.xml"), "709")
    smodel, hmodel = make_singlecrystal(verbose=True)

    Ts = np.array([500, 550, 600, 650.0]) + 273.15
    times = np.insert(np.logspace(0, 7, 100) * 3600.0, 0, [0])
    times = np.insert(times, 1, [1.0e2])

    for T in Ts:
        driver = drivers.Driver_sd(smodel, verbose=True, T_init=T)

        for t in times[1:]:
            print("t = {:.3E}".format(t))
            driver.strain_step(np.zeros((6,)), t, T)

        h = np.array(driver.stored_int)
        carbide1 = h[:, -6:-4]
        carbide2 = h[:, -4:-2]
        laves = h[:, -2:]

        hours = times / 3600.0

        r_carbide1 = carbide1[:, 0] * rs
        N_carbide1 = carbide1[:, 1] * Ns
        f_carbide1 = 4.0 / 3.0 * np.pi * N_carbide1 * r_carbide1**3
        A_carbide1 = 2.0 * r_carbide1 * N_carbide1

        r_carbide2 = carbide2[:, 0] * rs
        N_carbide2 = carbide2[:, 1] * Ns
        f_carbide2 = 4.0 / 3.0 * np.pi * N_carbide2 * r_carbide2**3
        A_carbide2 = 2.0 * r_carbide2 * N_carbide2

        r_laves = laves[:, 0] * rs
        N_laves = laves[:, 1] * Ns
        f_laves = 4.0 / 3.0 * np.pi * N_laves * r_laves**3
        A_laves = 2.0 * r_laves * N_laves

        taus = []
        for hi in h:
            blank = history.History()
            hist = history.History()
            smodel.populate_hist(hist)
            hist.set_data(hi)
            taus.append(hmodel.hist_to_tau(0, 0, hist, L, T, blank))

        plt.figure()
        plt.semilogx(hours, r_carbide1 * 1.0e9, label="Carbide 1")
        plt.semilogx(hours, r_carbide2 * 1.0e9, label="Carbide 2")
        plt.semilogx(hours, r_laves * 1.0e9, label="Laves")
        plt.legend(loc="best")
        plt.xlabel("Time (hours)")
        plt.ylabel("Radius (nm)")
        plt.xlim([1e0, 1e7])
        plt.show()

        plt.figure()
        plt.semilogx(hours, f_carbide1, label="Carbide 1")
        plt.semilogx(hours, f_carbide2, label="Carbide 2")
        plt.semilogx(hours, f_laves, label="Laves")
        plt.legend(loc="best")
        plt.xlabel("Time (hours)")
        plt.ylabel("Volume fraction")
        plt.xlim([1e0, 1e7])
        plt.show()

        plt.figure()
        plt.semilogx(hours, A_carbide1, label="Carbide 1")
        plt.semilogx(hours, A_carbide2, label="Carbide 2")
        plt.semilogx(hours, A_laves, label="Laves")
        plt.legend(loc="best")
        plt.xlabel("Time (hours)")
        plt.ylabel("Area density (1/m^2)")
        plt.xlim([1e0, 1e7])
        plt.show()

        plt.figure()
        plt.semilogx(hours, taus)
        plt.xlabel("Time (hours)")
        plt.ylabel("Slip resistance (MPa)")
        plt.xlim([1e0, 1e7])
        plt.show()
