from neml import models, interpolate, elasticity, history
from neml.cp import generalizedhucocks, hucocks, crystallography, sliprules, slipharden
from neml.math import rotations, tensors

import unittest
import numpy as np
import scipy.interpolate as inter

from common import differentiate_new
from nicediff import *
from test_slipharden import CommonSlipHardening


class TestPrecipitationModel(unittest.TestCase):
    def setUp(self):
        Ts = np.array([500.0, 550.0, 600.0, 650.0]) + 273.15

        self.L = crystallography.CubicLattice(1.0)
        self.L.add_slip_system([1, 1, 0], [1, 1, 1])

        self.Q = rotations.Orientation(35.0, 17.0, 14.0, angle_type="degrees")
        self.S = tensors.Symmetric(
            np.array([[100.0, -25.0, 10.0], [-25.0, -17.0, 15.0], [10.0, 15.0, 35.0]])
        )

        self.T = 600 + 273.15

        self.uf = 1.0

        self.J1 = 2e14 * self.uf**2.0
        self.J2 = 3.3e14 * self.uf**2.0
        self.K = 2.56e-30 / self.uf**4.0
        self.L0 = 3.1623e-7 / self.uf
        self.b = 2.5e-10
        self.bd = self.b / self.uf
        self.ad = 0.35
        self.G = interpolate.PiecewiseLinearInterpolate(
            list(Ts), [61068, 59541.0, 57633.6, 55725.2]
        )

        self.dmodel = hucocks.DislocationSpacingHardening(
            self.J1, self.J2, self.K, self.L0, self.ad, self.bd, self.G, self.L
        )

        # Chemical species
        # Suppose for 709 we need to track Cr, C, Mo, and Mn.
        # We assume there is ample Fe so we don't track it.
        self.Cr = generalizedhucocks.GeneralizedHuCocksSpecies(
            "Cr",
            0.1625,
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [0.1564, 0.1569, 0.1575, 0.1583]
            ),
        )
        self.C = generalizedhucocks.GeneralizedHuCocksSpecies(
            "C",
            0.000375,
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [7.25e-8, 2.92e-7, 9.48e-7, 2.97e-6]
            ),
        )
        self.Mo = generalizedhucocks.GeneralizedHuCocksSpecies(
            "Mo",
            0.0233,
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [0.0025, 0.0046, 0.0076, 0.0116]
            ),
        )
        self.Mn = generalizedhucocks.GeneralizedHuCocksSpecies(
            "Mn",
            0.1,
            0.05,
        )

        # Precipitates
        # Carbide1 Cr23C6
        # Carbide2 Mn3C
        # Laves Fe2Mo
        self.carbide1 = generalizedhucocks.GeneralizedHuCocksPrecipitate(
            "Cr23C6",  # composition
            "Cr C",  # species
            "Cr",  # the rate-limiting species
            [
                interpolate.PiecewiseLinearInterpolate(
                    list(Ts), [0.69845, 0.6905, 0.6832, 0.6752]
                ),
                interpolate.PiecewiseLinearInterpolate(
                    list(Ts), [0.0513, 0.0513, 0.0513, 0.0513]
                ),
            ],  # cp
            3.6e-10,  # am
            6e-6,  # Vm
            1.5e-4,  # D0
            240e3,  # Q0
            1e13,  # N0
            0.3,  # chi
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [1.0, 1.0, 0.3, 0.03]
            ),  # Cf
        )
        self.carbide2 = generalizedhucocks.GeneralizedHuCocksPrecipitate(
            "Mn3C",  # composition
            "Mn C",  # species
            "Mn",  # the rate-limiting species
            [
                interpolate.ConstantInterpolate(0.5),
                interpolate.ConstantInterpolate(0.1),
            ],  # cp
            3.6e-10,  # am
            5e-6,  # Vm
            1e-4,  # D0
            250e3,  # Q0
            5e13,  # N0
            0.25,  # chi
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [1.0, 1.0, 0.3, 0.03]
            ),  # Cf
        )
        self.laves = generalizedhucocks.GeneralizedHuCocksPrecipitate(
            "Fe2Mo",  # composition
            "Mo",  # species
            "Mo",  # the rate-limiting species
            [interpolate.ConstantInterpolate(0.5)],  # cp
            3.6e-10,  # am
            2e-6,  # Vm
            7.4e-4,  # D0
            283e3,  # Q0
            5e14,  # N0
            0.25,  # chi
            interpolate.ConstantInterpolate(1.0),  # Cf
        )

        # The precipitation model for 709 tracks the nucleation and growth of
        # 3 precipitates from 4 chemical species
        # Cr, C -> Cr23C6 (carbide1)
        # Mn, C -> Mn3C (carbide2)
        # Mo -> Fe2Mo (laves)
        self.pmodel = generalizedhucocks.GeneralizedHuCocksPrecipitationModel(
            [self.Cr, self.C, self.Mo, self.Mn],
            [self.carbide1, self.carbide2, self.laves],
        )

    def vector_state(self, vec):
        hist = history.History()

        hist.add_scalar("Cr23C6_r")
        hist.set_scalar("Cr23C6_r", vec[0])
        hist.add_scalar("Cr23C6_N")
        hist.set_scalar("Cr23C6_N", vec[1])

        hist.add_scalar("Mn3C_r")
        hist.set_scalar("Mn3C_r", vec[2])
        hist.add_scalar("Mn3C_N")
        hist.set_scalar("Mn3C_N", vec[3])

        hist.add_scalar("Fe2Mo_r")
        hist.set_scalar("Fe2Mo_r", vec[4])
        hist.add_scalar("Fe2Mo_N")
        hist.set_scalar("Fe2Mo_N", vec[5])

        return hist

    def test_history(self):
        """
        Paranoid check on the history state, given this is a kinda complicated model
        """
        H1 = history.History()
        self.pmodel.populate_hist(H1)
        self.pmodel.init_hist(H1)

        self.assertEqual(len(np.array(H1)), 3 + 3)

        should = np.array(
            [
                1.0e-9 / self.carbide1.rs(),
                1.0e11 / self.carbide1.Ns(),
            ]
            + [
                1.0e-9 / self.carbide2.rs(),
                1.0e11 / self.carbide2.Ns(),
            ]
            + [
                1.0e-9 / self.laves.rs(),
                1.0e11 / self.laves.Ns(),
            ]
        )

        self.assertTrue(np.allclose(np.array(H1), should))

    def test_jac_growth(self):
        """
        Test jacobian in the nucleation/growth regime
        """
        T = 550.0 + 273.15
        r = 1.0e-9
        N = 1.0e11

        state = np.array(
            [
                r / self.carbide1.rs(),
                N / self.carbide1.Ns(),
                r / self.carbide2.rs(),
                N / self.carbide2.Ns(),
                r / self.laves.rs(),
                N / self.laves.Ns(),
            ]
        )

        exact = self.pmodel.d_rate(T, self.vector_state(state))
        num = differentiate_new(
            lambda x: np.array(self.pmodel.rate(T, self.vector_state(x))),
            state,
        )

        self.assertTrue(np.allclose(exact, num))

    def test_jac_ripening(self):
        """
        Test jacobian in the ripening regime
        """
        T = 500.0 + 273.15
        r1 = 1.85864935e-07
        r2 = 9.19867846e-08
        r3 = 1.13153533e-07
        N1 = 1.05558617e17
        N2 = 7.19986332e17
        N3 = 4.98149587e18

        state = np.array(
            [
                r1 / self.carbide1.rs(),
                N1 / self.carbide1.Ns(),
                r2 / self.carbide2.rs(),
                N2 / self.carbide2.Ns(),
                r3 / self.laves.rs(),
                N3 / self.laves.Ns(),
            ]
        )

        exact = self.pmodel.d_rate(T, self.vector_state(state))
        num = differentiate_new(
            lambda x: np.array(self.pmodel.rate(T, self.vector_state(x))),
            state,
        )

        self.assertTrue(np.allclose(exact, num))


class TestGeneralizedHuCocksHardening(unittest.TestCase, CommonSlipHardening):
    def setUp(self):
        Ts = np.array([500.0, 550.0, 600.0, 650.0]) + 273.15

        self.L = crystallography.CubicLattice(1.0)
        self.L.add_slip_system([1, 1, 0], [1, 1, 1])

        self.Q = rotations.Orientation(35.0, 17.0, 14.0, angle_type="degrees")
        self.S = tensors.Symmetric(
            np.array([[100.0, -25.0, 10.0], [-25.0, -17.0, 15.0], [10.0, 15.0, 35.0]])
        )

        self.T = 600 + 273.15

        self.uf = 1.0

        self.J1 = 2e14 * self.uf**2.0
        self.J2 = 3.3e14 * self.uf**2.0
        self.K = 2.56e-30 / self.uf**4.0
        self.L0 = 3.1623e-7 / self.uf
        self.b = 2.5e-10
        self.bd = self.b / self.uf
        self.ad = 0.35
        self.G = interpolate.PiecewiseLinearInterpolate(
            list(Ts), [61068, 59541.0, 57633.6, 55725.2]
        )

        self.dmodel = hucocks.DislocationSpacingHardening(
            self.J1, self.J2, self.K, self.L0, self.ad, self.bd, self.G, self.L
        )

        # Chemical species
        # Suppose for 709 we need to track Cr, C, Mo, and Mn.
        # We assume there is ample Fe so we don't track it.
        self.Cr = generalizedhucocks.GeneralizedHuCocksSpecies(
            "Cr",
            0.1625,
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [0.1564, 0.1569, 0.1575, 0.1583]
            ),
        )
        self.C = generalizedhucocks.GeneralizedHuCocksSpecies(
            "C",
            0.000375,
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [7.25e-8, 2.92e-7, 9.48e-7, 2.97e-6]
            ),
        )
        self.Mo = generalizedhucocks.GeneralizedHuCocksSpecies(
            "Mo",
            0.0233,
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [0.0025, 0.0046, 0.0076, 0.0116]
            ),
        )
        self.Mn = generalizedhucocks.GeneralizedHuCocksSpecies(
            "Mn",
            0.1,
            0.05,
        )

        # Precipitates
        # Carbide1 Cr23C6
        # Carbide2 Mn3C
        # Laves Fe2Mo
        self.carbide1 = generalizedhucocks.GeneralizedHuCocksPrecipitate(
            "Cr23C6",  # composition
            "Cr C",  # species
            "Cr",  # the rate-limiting species
            [
                interpolate.PiecewiseLinearInterpolate(
                    list(Ts), [0.69845, 0.6905, 0.6832, 0.6752]
                ),
                interpolate.PiecewiseLinearInterpolate(
                    list(Ts), [0.0513, 0.0513, 0.0513, 0.0513]
                ),
            ],  # cp
            3.6e-10,  # am
            6e-6,  # Vm
            1.5e-4,  # D0
            240e3,  # Q0
            1e13,  # N0
            0.3,  # chi
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [1.0, 1.0, 0.3, 0.03]
            ),  # Cf
        )
        self.carbide2 = generalizedhucocks.GeneralizedHuCocksPrecipitate(
            "Mn3C",  # composition
            "Mn C",  # species
            "Mn",  # the rate-limiting species
            [
                interpolate.ConstantInterpolate(0.5),
                interpolate.ConstantInterpolate(0.5),
            ],  # cp
            3.6e-10,  # am
            5e-6,  # Vm
            1e-4,  # D0
            250e3,  # Q0
            5e13,  # N0
            0.25,  # chi
            interpolate.PiecewiseLinearInterpolate(
                list(Ts), [1.0, 1.0, 0.3, 0.03]
            ),  # Cf
        )
        self.laves = generalizedhucocks.GeneralizedHuCocksPrecipitate(
            "Fe2Mo",  # composition
            "Mo",  # species
            "Mo",  # the rate-limiting species
            [interpolate.ConstantInterpolate(0.5)],  # cp
            3.6e-10,  # am
            2e-6,  # Vm
            7.4e-4,  # D0
            283e3,  # Q0
            5e14,  # N0
            0.25,  # chi
            interpolate.ConstantInterpolate(1.0),  # Cf
        )

        # The precipitation model for 709 tracks the nucleation and growth of
        # 3 precipitates from 4 chemical species
        # Cr, C -> Cr23C6 (carbide1)
        # Mn, C -> Mn3C (carbide2)
        # Mo -> Fe2Mo (laves)
        self.pmodel = generalizedhucocks.GeneralizedHuCocksPrecipitationModel(
            [self.Cr, self.C, self.Mo, self.Mn],
            [self.carbide1, self.carbide2, self.laves],
        )

        self.ap = 0.84
        self.ac = 0.000457
        self.b = 2.5e-10
        self.G = interpolate.PiecewiseLinearInterpolate(
            list(Ts), [61068, 59541.0, 57633.6, 55725.2]
        )
        self.model = generalizedhucocks.GeneralizedHuCocksHardening(
            self.dmodel, self.pmodel, self.ap, self.ac, self.b, self.G
        )

        self.g0 = 1.0
        self.n = 3.0
        self.sliprule = sliprules.PowerLawSlipRule(self.model, self.g0, self.n)
        self.nslip = self.L.ntotal

        self.current = 1.15e-7
        self.r_carbide1 = 1e-9 / self.carbide1.rs()
        self.N_carbide1 = 1e11 / self.carbide1.Ns()
        self.r_carbide2 = 1e-9 / self.carbide2.rs()
        self.N_carbide2 = 1e11 / self.carbide2.Ns()
        self.r_laves = 1e-9 / self.laves.rs()
        self.N_laves = 1e11 / self.laves.Ns()

        self.H = history.History()
        for i in range(self.nslip):
            self.H.add_scalar("spacing_" + str(i))
            self.H.set_scalar("spacing_" + str(i), self.current)

        self.H.add_scalar("Cr23C6_r")
        self.H.set_scalar("Cr23C6_r", self.r_carbide1)
        self.H.add_scalar("Cr23C6_N")
        self.H.set_scalar("Cr23C6_N", self.N_carbide1)

        self.H.add_scalar("Mn3C_r")
        self.H.set_scalar("Mn3C_r", self.r_carbide2)
        self.H.add_scalar("Mn3C_N")
        self.H.set_scalar("Mn3C_N", self.N_carbide2)

        self.H.add_scalar("Fe2Mo_r")
        self.H.set_scalar("Fe2Mo_r", self.r_laves)
        self.H.add_scalar("Fe2Mo_N")
        self.H.set_scalar("Fe2Mo_N", self.N_laves)

        self.fixed = history.History()

    def test_history(self):
        """
        Paranoid check on the history state, given this is a kinda complicated model
        """
        H1 = history.History()
        self.model.populate_hist(H1)
        self.model.init_hist(H1)

        self.assertEqual(len(np.array(H1)), 12 + 3 + 3)

        should = np.array(
            [self.L0] * 12
            + [
                1.0e-9 / self.carbide1.rs(),
                1.0e11 / self.carbide1.Ns(),
            ]
            + [
                1.0e-9 / self.carbide2.rs(),
                1.0e11 / self.carbide2.Ns(),
            ]
            + [
                1.0e-9 / self.laves.rs(),
                1.0e11 / self.laves.Ns(),
            ]
        )

        self.assertTrue(np.allclose(np.array(H1), should))
