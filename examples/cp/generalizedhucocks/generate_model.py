import numpy as np
from neml import interpolate, elasticity, history
from neml.cp import (
    hucocks,
    generalizedhucocks,
    crystallography,
    inelasticity,
    kinematics,
    singlecrystal,
)
from neml.math import rotations
import sys

sys.path.append("../../..")


def make_pmodel():
    Ts = np.array([500.0, 550.0, 600.0, 650.0]) + 273.15

    # Chemical species
    # Suppose for 709 we need to track Cr, C, Mo, and Mn.
    # We assume there is ample Fe so we don't track it.
    Cr = generalizedhucocks.GeneralizedHuCocksSpecies(
        "Cr",
        0.1625,
        interpolate.PiecewiseLinearInterpolate(
            list(Ts), [0.1564, 0.1569, 0.1575, 0.1583]
        ),
    )
    C = generalizedhucocks.GeneralizedHuCocksSpecies(
        "C",
        0.000375,
        interpolate.PiecewiseLinearInterpolate(
            list(Ts), [7.25e-8, 2.92e-7, 9.48e-7, 2.97e-6]
        ),
    )
    Mo = generalizedhucocks.GeneralizedHuCocksSpecies(
        "Mo",
        0.0233,
        interpolate.PiecewiseLinearInterpolate(
            list(Ts), [0.0025, 0.0046, 0.0076, 0.0116]
        ),
    )
    Mn = generalizedhucocks.GeneralizedHuCocksSpecies(
        "Mn",
        0.1,
        0.05,
    )

    # Precipitates
    # Carbide1 Cr23C6
    # Carbide2 Mn3C
    # Laves Fe2Mo
    carbide1 = generalizedhucocks.GeneralizedHuCocksPrecipitate(
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
        interpolate.PiecewiseLinearInterpolate(list(Ts), [1.0, 1.0, 0.3, 0.03]),  # Cf
    )
    carbide2 = generalizedhucocks.GeneralizedHuCocksPrecipitate(
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
        interpolate.PiecewiseLinearInterpolate(list(Ts), [1.0, 1.0, 0.3, 0.03]),  # Cf
    )
    laves = generalizedhucocks.GeneralizedHuCocksPrecipitate(
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
    pmodel = generalizedhucocks.GeneralizedHuCocksPrecipitationModel(
        [Cr, C, Mo, Mn], [carbide1, carbide2, laves]
    )

    return pmodel


def make_singlecrystal(verbose=False):
    Ts = np.array([500.0, 550.0, 600.0, 650.0]) + 273.15

    L = crystallography.CubicLattice(1.0)
    L.add_slip_system([1, 1, 0], [1, 1, 1])

    uf = 1.0e-9

    J1 = 2e14 * uf**2.0
    J2 = 3.3e14 * uf**2.0
    K = 2.56e-30 / uf**4.0
    L0 = 3.1623e-7 / uf
    b = 2.5e-10
    b_d = b / uf
    ad = 0.35
    G = interpolate.PiecewiseLinearInterpolate(
        list(Ts), [61068, 59541.0, 57633.6, 55725.2]
    )

    dmodel = hucocks.DislocationSpacingHardening(J1, J2, K, L0, ad, b_d, G, L)

    pmodel = make_pmodel()

    ap = 0.84
    ac = 0.000457

    tau_model = generalizedhucocks.GeneralizedHuCocksHardening(
        dmodel, pmodel, ap, ac, b, G
    )

    g0 = 1.0
    a0 = 0.5
    G0 = 77000.0e6
    A = 3.0 / 4.0
    B = 4.0 / 3.0

    slip_model = hucocks.ArrheniusSlipRule(tau_model, g0, A, B, b, a0, G0)

    imodel = inelasticity.AsaroInelasticity(slip_model)

    youngs = interpolate.PiecewiseLinearInterpolate(
        list(Ts), [160000.0, 156000.0, 151000.0, 140000.0]
    )
    emodel = elasticity.IsotropicLinearElasticModel(youngs, "youngs", 0.31, "poissons")

    kmodel = kinematics.StandardKinematicModel(emodel, imodel)

    smodel = singlecrystal.SingleCrystalModel(
        kmodel,
        L,
        linesearch=True,
        initial_rotation=rotations.CrystalOrientation(0, 0, 0, angle_type="degrees"),
        verbose=verbose,
    )

    return smodel, tau_model
