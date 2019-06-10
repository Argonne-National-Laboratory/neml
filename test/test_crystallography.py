#!/usr/bin/env python

from neml.math import tensors, rotations
from neml.cp import crystallography

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
h = -1.0/2.0
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

class TestRotationalGroups(unittest.TestCase):
  def compare(self, matrix, quaternion):
    for M, q in zip(matrix, quaternion):
      self.assertTrue(np.allclose(M, q.to_matrix()))

  def test_1(self):
    self.compare(ktw_rotational_groups("1"), 
        crystallography.symmetry_rotations("1"))

  def test_2(self):
    self.compare(ktw_rotational_groups("2"),
        crystallography.symmetry_rotations("2"))

  def test_222(self):
    self.compare(ktw_rotational_groups("222"),
        crystallography.symmetry_rotations("2"))

  def test_42(self):
    self.compare(ktw_rotational_groups("42"),
        crystallography.symmetry_rotations("42"))

  def test_4(self):
    self.compare(ktw_rotational_groups("4"),
        crystallography.symmetry_rotations("4"))

  def test_3(self):
    self.compare(ktw_rotational_groups("3"),
        crystallography.symmetry_rotations("3"))

  def test_6(self):
    self.compare(ktw_rotational_groups("6"),
        crystallography.symmetry_rotations("6"))

  def test_32(self):
    self.compare(ktw_rotational_groups("32"),
        crystallography.symmetry_rotations("32"))

  def test_622(self):
    self.compare(ktw_rotational_groups("622"),
        crystallography.symmetry_rotations("622"))

  def test_23(self):
    self.compare(ktw_rotational_groups("23"),
        crystallography.symmetry_rotations("23"))

  def test_432(self):
    self.compare(ktw_rotational_groups("432"),
        crystallography.symmetry_rotations("432"))

class LTests(object):
  def test_lattice(self):
    test = np.zeros((3,3))

    avs = [self.lattice.a1, self.lattice.a2, self.lattice.a3]
    bvs = [self.lattice.b1, self.lattice.b2, self.lattice.b3]

    for i in range(3):
      for j in range(3):
        test[i,j] = avs[i].dot(bvs[j])
    
    self.assertTrue(np.allclose(np.eye(3), test))

class TestCubic(unittest.TestCase, LTests):
  def setUp(self):
    self.a = 1.3
    self.lattice = crystallography.CubicLattice(self.a)

  def test_fcc(self):
    self.lattice.add_slip_system([1,1,0],[1,1,1])

    self.assertEqual(len(self.lattice.burgers_vectors[0]), 12)
    self.assertEqual(len(self.lattice.slip_directions[0]), 12)
    self.assertEqual(len(self.lattice.slip_planes[0]), 12)

    for d,n in zip(self.lattice.slip_directions[0], self.lattice.slip_planes[0]):
      self.assertAlmostEqual(d.norm(), 1.0)
      self.assertAlmostEqual(n.norm(), 1.0)
      self.assertAlmostEqual(d.dot(n), 0.0)
    
    for b in self.lattice.burgers_vectors[0]:
      self.assertAlmostEqual(b.norm(), np.sqrt(2.0) * self.a)

  def test_bcc(self):
    self.lattice.add_slip_system([1,1,1],[1,1,0])
    self.lattice.add_slip_system([1,1,1],[1,2,3])
    self.lattice.add_slip_system([1,1,1],[1,1,2])
    
    self.assertEqual(len(self.lattice.slip_directions[0]), 12)
    self.assertEqual(len(self.lattice.slip_directions[1]), 24)
    self.assertEqual(len(self.lattice.slip_directions[2]), 12)

    for dg,ng,bg in zip(self.lattice.slip_directions, self.lattice.slip_planes,
        self.lattice.burgers_vectors):
      for d,n in zip(dg, ng):
        self.assertAlmostEqual(d.norm(), 1.0)
        self.assertAlmostEqual(n.norm(), 1.0)
        self.assertAlmostEqual(d.dot(n), 0.0)

      for b in bg:
        self.assertAlmostEqual(b.norm(), np.sqrt(3.0) * self.a)
