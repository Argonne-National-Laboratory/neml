#!/usr/bin/env python3

import unittest

import numpy as np

from neml import history, elasticity
from neml.math import tensors, rotations, matrix
from neml.cp import crystallography, slipharden, sliprules, postprocessors, singlecrystal, inelasticity, kinematics

from common import *

class TestPTRTwinReorientation(unittest.TestCase):
  def setUp(self):
    a = 2.9511*0.1 # nm
    c = 4.68433*0.1 # nm

    C11 = 160000.0
    C33 = 181000.0
    C44 = 46500.0
    C12 = 90000.0
    C13 = 66000.0

    tau0 = np.array([170.0]*3+[90.5]*3+[210]*6+[180.0]*6+[250.0]*6)
    
    # Not realistic but I just want it to roll quickly
    H1 = 100.0 
    H2 = 100.0

    g0 = 1.0
    n = 12.0

    M = matrix.SquareMatrix(24, type = "diagonal_blocks", 
        data = [H1,H2], blocks = [12,12])

    self.lattice = crystallography.HCPLattice(a, c)
    # Basal <a>
    self.lattice.add_slip_system([1,1,-2,0],[0,0,0,1])
    # Prismatic <a>
    self.lattice.add_slip_system([1,1,-2,0],[1,0,-1,0])
    # Pyramidal <c+a>
    self.lattice.add_slip_system([1,1,-2,3],[1,1,-2,2])
    # Tension twinning
    self.lattice.add_twin_system([-1,0,1,1],[1,0,-1,2],[1,0,-1,1],[1,0,-1,-2])
    # Compression twinning
    self.lattice.add_twin_system([1,1,-2,-3],[1,1,-2,2],[2,2,-4,3],[1,1,-2,-4])

    emodel = elasticity.TransverseIsotropicLinearElasticModel(
        C11,C33,C12,C13,C44,"components")

    strength = slipharden.SimpleLinearHardening(M, tau0)
    slipmodel = sliprules.PowerLawSlipRule(strength, g0, n)
    imodel = inelasticity.AsaroInelasticity(slipmodel)
    kmodel = kinematics.StandardKinematicModel(emodel, imodel)

    self.Q = rotations.CrystalOrientation(35.0,17.0,14.0, angle_type = "degrees")
    
    L = np.array([[0.5,0.5,0],[0,0.5,0],[0,0,-1.0]])
    self.D = tensors.Symmetric(0.5*(L+L.T))
    self.W = tensors.Skew(0.5*(L-L.T))

    self.model = singlecrystal.SingleCrystalModel(kmodel, self.lattice,
        initial_rotation = self.Q)

    self.threshold = 0.1
    self.processor = postprocessors.PTRTwinReorientation(self.threshold)

  def test_history(self):
    h = history.History()
    self.processor.populate_history(self.lattice, h)
    self.processor.init_history(self.lattice, h)

    names = ["twin_fraction"+str(i) for i in range(12,24)] + ["twinned"]

    self.assertEqual(names, h.items)

    self.assertTrue(np.allclose(h, np.zeros((13,))))
   

