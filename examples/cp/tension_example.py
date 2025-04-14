#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np
import numpy.random as ra

from neml.cp import polycrystal, crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polefigures
from neml.math import rotations, tensors, nemlmath
from neml import elasticity, drivers

import matplotlib.pyplot as plt

if __name__ == "__main__":
    nthreads = 32
    ngrain = 500

    # Let's pretend these are active convention
    orientations = rotations.random_orientations(ngrain) 

    # For plotting get in the passive convention
    orientations_passive = [o.inverse() for o in orientations]

    # Setup the model with fairly arbitrary properties
    C11 = 287750.0  # MPa
    C12 = 127920.0   # MPa
    C44 = 120710.0   # MPa

    t0 = 50.0
    ts = 40.0
    b = 1.0

    g0 = 1.0
    n = 12.0

    lattice = crystallography.CubicLattice(1.0)
    lattice.add_slip_system([1,1,0],[1,1,1])

    strengthmodel = slipharden.VoceSlipHardening(ts, b, t0)
    slipmodel = sliprules.PowerLawSlipRule(strengthmodel, g0, n)
    imodel = inelasticity.AsaroInelasticity(slipmodel)
    emodel = elasticity.CubicLinearElasticModel(C11,C12,C44, "components")
    kmodel = kinematics.StandardKinematicModel(emodel, imodel)

    model = singlecrystal.SingleCrystalModel(kmodel, lattice)

    pmodel = polycrystal.TaylorModel(model, orientations, nthreads = nthreads)

   
    # Pass in passive convention, ugh, sorry
    polefigures.inverse_pole_figure_discrete(orientations_passive, [1,0,0], lattice, reduce_figure = "cubic", axis_labels = ["100", "110", "111"])
    plt.show()

    # Strain rate and target strain for loading
    erate = 8.33e-5
    emax = 0.01
    nsteps = 50

    res = drivers.uniaxial_test(pmodel, erate, emax = emax, nsteps = nsteps, 
                                full_results = True, verbose = True)

    internal_state = np.array(res['history'])
    
    # Final orientations, in the passive convention
    final_orientations = pmodel.orientations(internal_state[-1])
    
    # Stress history for each crystal
    nhist = model.nstore
    stress = internal_state[:,pmodel.n*nhist:pmodel.n*nhist+6*pmodel.n:6]

    # Overall strain history
    strain = np.array(res['strain'])[:,0]

    # Overall stress
    stress_overall = np.array(res['stress'])[:,0]
    
    # Pass in passive convention, ugh, sorry
    polefigures.inverse_pole_figure_discrete(final_orientations, [1,0,0], lattice, reduce_figure = "cubic", axis_labels = ["100", "110", "111"],
            color = stress[-1])
    plt.show()
