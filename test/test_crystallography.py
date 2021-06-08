#!/usr/bin/env python

from neml.math import tensors, rotations
from neml.cp import crystallography

from common import differentiate

import unittest
import numpy as np
import numpy.linalg as la

ktw_tetra = np.array(
    [
    [ [ 1, 0, 0],
      [ 0, 1, 0],
      [ 0, 0, 1] ],
    [ [-1, 0, 0],
      [ 0, 1, 0],
      [ 0, 0,-1] ],
    [ [ 1, 0, 0],
      [ 0,-1, 0],
      [ 0, 0,-1] ],
    [ [-1, 0, 0],
      [ 0,-1, 0],
      [ 0, 0, 1] ],
    [ [ 0, 1, 0],
      [-1, 0, 0],
      [ 0, 0, 1] ],
    [ [ 0,-1, 0],
      [ 1, 0, 0],
      [ 0, 0, 1] ],
    [ [ 0, 1, 0],
      [ 1, 0, 0],
      [ 0, 0,-1] ],
    [ [ 0,-1, 0],
      [-1, 0, 0],
      [ 0, 0,-1] ]
    ], dtype = float)

a = np.sqrt(3.0)/2.0
h = 1.0/2.0
ktw_hexagonal = np.array(
    [
      [ [ 1, 0, 0], 
        [ 0, 1, 0],
        [ 0, 0, 1] ],
      [ [-h, a, 0],
        [-a,-h, 0],
        [ 0, 0, 1] ],
      [ [-h,-a, 0],
        [ a,-h, 0],
        [ 0, 0, 1] ],
      [ [ h, a, 0],
        [-a, h, 0],
        [ 0, 0, 1] ],
      [ [-1, 0, 0],
        [ 0,-1, 0],
        [ 0, 0, 1] ],
      [ [ h,-a, 0],
        [ a, h, 0],
        [ 0, 0, 1] ],
      [ [-h,-a, 0],
        [-a, h, 0],
        [ 0, 0,-1] ],
      [ [ 1, 0, 0],
        [ 0,-1, 0],
        [ 0, 0,-1] ],
      [ [-h, a, 0],
        [ a, h, 0],
        [ 0, 0,-1] ],
      [ [ h, a, 0],
        [ a,-h, 0],
        [ 0, 0,-1] ],
      [ [-1, 0, 0],
        [ 0, 1, 0],
        [ 0, 0,-1] ],
      [ [ h,-a, 0],
        [-a,-h, 0],
        [ 0, 0,-1 ] ]
      ])

ktw_cubic = np.array(
    [
      [ [ 1, 0, 0],
        [ 0, 1, 0],
        [ 0, 0, 1] ],
      [ [ 0, 0, 1],
        [ 1, 0, 0],
        [ 0, 1, 0] ],
      [ [ 0, 1, 0],
        [ 0, 0, 1],
        [ 1, 0, 0] ],
      [ [ 0,-1, 0],
        [ 0, 0, 1],
        [-1, 0, 0] ],
      [ [ 0,-1, 0],
        [ 0, 0,-1],
        [ 1, 0, 0] ],
      [ [ 0, 1, 0],
        [ 0, 0,-1],
        [-1, 0, 0] ],
      [ [ 0, 0,-1],
        [ 1, 0, 0],
        [ 0,-1, 0] ],
      [ [ 0, 0,-1],
        [-1, 0, 0],
        [ 0, 1, 0] ],
      [ [ 0, 0, 1],
        [-1, 0, 0],
        [ 0,-1, 0] ],
      [ [-1, 0, 0],
        [ 0, 1, 0],
        [ 0, 0,-1] ],
      [ [-1, 0, 0],
        [ 0,-1, 0],
        [ 0, 0, 1] ],
      [ [ 1, 0, 0],
        [ 0,-1, 0],
        [ 0, 0,-1] ],
      [ [ 0, 0,-1],
        [ 0,-1, 0],
        [-1, 0, 0] ],
      [ [ 0, 0, 1],
        [ 0,-1, 0],
        [ 1, 0, 0] ],
      [ [ 0, 0, 1],
        [ 0, 1, 0],
        [-1, 0, 0] ],
      [ [ 0, 0,-1],
        [ 0, 1, 0],
        [ 1, 0, 0] ],
      [ [-1, 0, 0],
        [ 0, 0,-1],
        [ 0,-1, 0] ],
      [ [ 1, 0, 0],
        [ 0, 0,-1],
        [ 0, 1, 0] ],
      [ [ 1, 0, 0],
        [ 0, 0, 1],
        [ 0,-1, 0] ],
      [ [-1, 0, 0],
        [ 0, 0, 1],
        [ 0, 1, 0] ],
      [ [ 0,-1, 0],
        [-1, 0, 0],
        [ 0, 0,-1] ],
      [ [ 0, 1, 0],
        [-1, 0, 0],
        [ 0, 0, 1] ],
      [ [ 0, 1, 0],
        [ 1, 0, 0],
        [ 0, 0,-1] ],
      [ [ 0,-1, 0],
        [ 1, 0, 0],
        [ 0, 0, 1] ]
      ])

def ktw_rotational_groups(typ):
  if typ == "1":
    return ktw_tetra[0:1]
  elif typ == "2":
    return ktw_tetra[0:2]
  elif typ == "222":
    return ktw_tetra[0:4]
  elif typ == "42":
    return ktw_tetra[:]
  elif typ == "4":
    return ktw_tetra[[0,3,4,5]]
  elif typ == "3":
    return ktw_hexagonal[0:3]
  elif typ == "6":
    return ktw_hexagonal[0:6]
  elif typ == "32":
    return ktw_hexagonal[[0,1,2,9,10,11]]
  elif typ == "622":
    return ktw_hexagonal[:]
  elif typ == "23":
    return ktw_cubic[0:12]
  elif typ == "432":
    return ktw_cubic[:]
  else:
    raise ValueError("Unknown rotational group %s." % typ)

groups = ["1", "2", "222", "42", "4", "3", "6", "32", "622", "23", "432"]

class TestRotationalGroups(unittest.TestCase):
  def compare(self, matrix, quaternion):
    for i,(M, q) in enumerate(zip(matrix, quaternion)):
      self.assertTrue(np.allclose(M, q.to_matrix()))

  def test_definition(self):
    for grp in groups:
      print(grp)
      self.compare(ktw_rotational_groups(grp),
          crystallography.symmetry_rotations(grp))

  def test_is_group(self):
    for grp in groups:
      print(grp)
      ops = crystallography.symmetry_rotations(grp)
      self.assertTrue(self.includes_identity(ops))
      self.assertTrue(self.closed(ops))
      self.assertTrue(self.has_inverse(ops))

  def includes_identity(self, ops):
    for op in ops:
      if np.allclose(op.quat, [1,0,0,0]):
        break
    else:
      return False
    return True

  def closed(self, ops):
    for i,op_a in enumerate(ops):
      for j,op_b in enumerate(ops):
        c = op_a * op_b
        found = False
        for op_c in ops:
          if np.allclose(c.quat, op_c.quat) or np.allclose(c.quat, -op_c.quat):
            found = True
            break
        if not found:
          return False

    return True

  def has_inverse(self, ops):
    for op_a in ops:
      opa_inv = op_a.inverse()
      for op_b in ops:
        if np.allclose(opa_inv.quat, op_b.quat) or np.allclose(opa_inv.quat, -op_b.quat):
          break
      else:
        return False
    
    return True

class LTests(object):
  def test_lattice(self):
    test = np.zeros((3,3))

    avs = [self.lattice.a1, self.lattice.a2, self.lattice.a3]
    bvs = [self.lattice.b1, self.lattice.b2, self.lattice.b3]

    for i in range(3):
      for j in range(3):
        test[i,j] = avs[i].dot(bvs[j])
    
    self.assertTrue(np.allclose(np.eye(3), test))

class ShearTests(object):
  def test_M(self):
    for i in range(self.lattice.ngroup):
      for j in range(self.lattice.nslip(i)):
        self.assertEqual(
            tensors.Symmetric(
              np.dot(self.QM,np.dot(np.outer(self.lattice.slip_directions[i][j].data,
                self.lattice.slip_planes[i][j].data), self.QM.T))),
              self.lattice.M(i,j,self.Q))

  def test_W(self):
    for i in range(self.lattice.ngroup):
      for j in range(self.lattice.nslip(i)):
        self.assertEqual(
            tensors.Skew(
              np.dot(self.QM,np.dot(np.outer(self.lattice.slip_directions[i][j].data,
                self.lattice.slip_planes[i][j].data), self.QM.T))),
              self.lattice.N(i,j,self.Q))

  def test_shear(self):
    for i in range(self.lattice.ngroup):
      for j in range(self.lattice.nslip(i)):
        self.assertTrue(np.isclose(
            np.einsum('ij,ij',
              np.dot(self.QM,np.dot(np.outer(self.lattice.slip_directions[i][j].data,
                self.lattice.slip_planes[i][j].data), self.QM.T)), self.S),
              self.lattice.shear(i,j,self.Q,self.ST)))

  def test_dshear(self):
    for i in range(self.lattice.ngroup):
      for j in range(self.lattice.nslip(i)):
        rs = lambda S: np.einsum('ij,ij',
              np.dot(self.QM,np.dot(np.outer(self.lattice.slip_directions[i][j].data,
                self.lattice.slip_planes[i][j].data), self.QM.T)), S)
        num = tensors.Symmetric(differentiate(rs, self.S)[0])

class TestCubicFCC(unittest.TestCase, LTests, ShearTests):
  def setUp(self):
    self.a = 1.3
    self.lattice = crystallography.CubicLattice(self.a)

    self.S = np.array([[20.0,-15.0,12.0],[-15.0,-40.0,5.0],[12.0,5.0,60.0]])
    self.ST = tensors.Symmetric(self.S)
    self.Q = rotations.Orientation(30.0,43.0,10.0, angle_type = "degrees")
    self.QM = self.Q.to_matrix()

    self.lattice.add_slip_system([1,1,0],[1,1,1])

  def test_fcc(self):
    self.assertEqual(len(self.lattice.burgers_vectors[0]), 12)
    self.assertEqual(len(self.lattice.slip_directions[0]), 12)
    self.assertEqual(len(self.lattice.slip_planes[0]), 12)

    self.assertEqual(self.lattice.ngroup, 1)
    self.assertEqual(self.lattice.nslip(0), 12)
    self.assertEqual(self.lattice.ntotal, 12)

    for d,n in zip(self.lattice.slip_directions[0], self.lattice.slip_planes[0]):
      self.assertAlmostEqual(d.norm(), 1.0)
      self.assertAlmostEqual(n.norm(), 1.0)
      self.assertAlmostEqual(d.dot(n), 0.0)
    
    for b in self.lattice.burgers_vectors[0]:
      self.assertAlmostEqual(b.norm(), np.sqrt(2.0) * self.a)

  def test_planes(self):
    self.assertEqual(self.lattice.nplanes, 4)
    for i in range(12):
      from_normals = self.lattice.unique_planes[self.lattice.plane_index(0, i)]
      from_direct = self.lattice.slip_planes[0][i]
      self.assertAlmostEqual(np.abs(from_normals.dot(from_direct)), 1.0)

class TestCubicBCC(unittest.TestCase, LTests, ShearTests):
  def setUp(self):
    self.a = 1.3
    self.lattice = crystallography.CubicLattice(self.a)

    self.S = np.array([[20.0,-15.0,12.0],[-15.0,-40.0,5.0],[12.0,5.0,60.0]])
    self.ST = tensors.Symmetric(self.S)
    self.Q = rotations.Orientation(30.0,43.0,10.0, angle_type = "degrees")
    self.QM = self.Q.to_matrix()

    self.lattice.add_slip_system([1,1,1],[1,1,0])
    self.lattice.add_slip_system([1,1,1],[1,2,3])
    self.lattice.add_slip_system([1,1,1],[1,1,2])

  def test_bcc(self):
    self.assertEqual(len(self.lattice.slip_directions[0]), 12)
    self.assertEqual(len(self.lattice.slip_directions[1]), 24)
    self.assertEqual(len(self.lattice.slip_directions[2]), 12)

    self.assertEqual(self.lattice.ngroup, 3)
    self.assertEqual(self.lattice.nslip(0), 12)
    self.assertEqual(self.lattice.nslip(1), 24)
    self.assertEqual(self.lattice.nslip(2), 12)
    self.assertEqual(self.lattice.ntotal, 48)

    for dg,ng,bg in zip(self.lattice.slip_directions, self.lattice.slip_planes,
        self.lattice.burgers_vectors):
      for d,n in zip(dg, ng):
        self.assertAlmostEqual(d.norm(), 1.0)
        self.assertAlmostEqual(n.norm(), 1.0)
        self.assertAlmostEqual(d.dot(n), 0.0)

      for b in bg:
        self.assertAlmostEqual(b.norm(), np.sqrt(3.0) * self.a)

  def test_planes(self):
    self.assertEqual(self.lattice.nplanes, 42)
    for g in range(self.lattice.ngroup):
      for i in range(self.lattice.nslip(g)):
        from_normals = self.lattice.unique_planes[self.lattice.plane_index(g, i)]
        from_direct = self.lattice.slip_planes[g][i]
        self.assertAlmostEqual(np.abs(from_normals.dot(from_direct)), 1.0)
