#!/usr/bin/env python

from neml import elasticity
from neml.math import tensors, projections

import common

import unittest
import numpy as np
import numpy.linalg as la

class TestProjectionsOnR2(unittest.TestCase):
  def setUp(self):
    self.s_full = tensors.RankTwo(np.array([[1,2,3],[4,5,6],[7,8,9.0]]))
    self.s_sym = tensors.Symmetric(np.array([[1,2,3],[2,4,5],[3,5,6]]))

  def test_full_normal(self):
    v = tensors.Vector([0,1.0,0])
    P = projections.normal_projection(v)
    
    proj = P.dot(self.s_full)

    self.assertEqual(proj, tensors.RankTwo(np.array([[0,0,0],[0,5,0],[0,0,0]])))

  def test_sym_normal(self):
    v = tensors.Vector([0,1.0,0])
    P = projections.normal_projection_ss(v)

    proj = P.dot(self.s_sym)

    self.assertEqual(proj, tensors.Symmetric(np.array([[0,0,0],[0,4,0],[0,0,0]])))
  
  def test_full_shear(self):
    v = tensors.Vector([0,1.0,0])
    P = projections.shear_projection(v)
    
    proj = P.dot(self.s_full)

    self.assertEqual(proj, tensors.RankTwo(np.array([[0,2,0],[0,0,0],[0,8,0]])))

  def test_sym_shear(self):
    v = tensors.Vector([0,1.0,0])
    P = projections.shear_projection_ss(v)

    proj = P.dot(self.s_sym)
    self.assertEqual(proj, tensors.Symmetric([[0,2,0],[2,0,5],[0,5,0]]))

  def test_normal_usecase(self):
    """
      Actually check that the stiffness in the appropriate direction resolves
      to zero
    """
    strain =  tensors.Symmetric(np.array([[1,1,1],[1,1,1],[1,1,1]]))
    emodel = elasticity.CubicLinearElasticModel(100000.0,0.3,60000.0, "moduli")
    C = emodel.C_tensor(0)
    
    n = tensors.Vector([1.0,-1,1]).normalize()

    P = projections.normal_projection_ss(n)
    PP = tensors.SymSymR4.id() - P
    PC = PP.dot(C)
    stress = PC.dot(strain).to_full()
    val = stress.dot(n).dot(n)

    self.assertAlmostEqual(val, 0)

    # Some random orthogonal direction
    s = n.cross(tensors.Vector([1,0.0,0])).normalize()
    val = stress.dot(s).dot(s)
    val2 = C.dot(strain).dot(s).dot(s)
    self.assertAlmostEqual(val,val2)


  def test_shear_usecase(self):
    """
      Actually check that the stiffness in the shear directions resolves to zero
    """
    strain =  tensors.Symmetric(np.array([[1,1,1],[1,1,1],[1,1,1]]))
    emodel = elasticity.CubicLinearElasticModel(100000.0,0.3,60000.0, "moduli")
    C = emodel.C_tensor(0)
    
    n = tensors.Vector([1.0,-1,1]).normalize()
    s = n.cross(tensors.Vector([1,0,0])).normalize()
    t = n.cross(s)

    P = projections.shear_projection_ss(n)
    PP = tensors.SymSymR4.id() - P
    PC = PP.dot(C)
    stress = PC.dot(strain).to_full()

    v1 = stress.dot(n).dot(s)
    self.assertAlmostEqual(v1,0)

    v2 = stress.dot(n).dot(t)
    self.assertAlmostEqual(v2,0)

    extra = (s+t).normalize()
    v3 = stress.dot(n).dot(extra)
    self.assertAlmostEqual(v3,0)
    
    # IDK, i guess check normal on a different plane
    vn1 = stress.dot(n).dot(n)
    vn2 = C.dot(strain).dot(n).dot(n)
    self.assertAlmostEqual(vn1,vn2)

    # IDK, i guess check shear on a different plane
    vn1 = stress.dot(s).dot(t)
    vn2 = C.dot(strain).dot(s).dot(t)
    self.assertAlmostEqual(vn1,vn2)
