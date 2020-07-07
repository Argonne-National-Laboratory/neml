#!/usr/bin/env python

from neml import history
from neml.math import rotations, tensors

import unittest

import numpy as np

class TestBasicAddGet(unittest.TestCase):
  def setUp(self):
    self.hist = history.History()
    self.hist.add_scalar("a")
    self.hist.add_vector("b")
    self.hist.add_ranktwo("C")
    self.hist.add_symmetric("D")
    self.hist.add_skew("E")
    self.hist.add_orientation("q")

  def test_exists(self):
    self.assertTrue(self.hist.contains("a"))
    self.assertFalse(self.hist.contains("alsd"))

  def test_get_correct(self):
    a = self.hist.get_scalar("a")
    b = self.hist.get_vector("b")
    C = self.hist.get_ranktwo("C")
    D = self.hist.get_symmetric("D")
    E = self.hist.get_skew("E")
    q = self.hist.get_orientation("q")

  def test_get_wrong_scalar(self):
    with self.assertRaises(RuntimeError):
      a = self.hist.get_scalar("b")

  def test_get_wrong_vector(self):
    with self.assertRaises(RuntimeError):
      a = self.hist.get_vector("a")

  def test_get_wrong_ranktwo(self):
    with self.assertRaises(RuntimeError):
      a = self.hist.get_ranktwo("b")

  def test_get_wrong_symmetric(self):
    with self.assertRaises(RuntimeError):
      a = self.hist.get_symmetric("b")

  def test_get_wrong_skew(self):
    with self.assertRaises(RuntimeError):
      a = self.hist.get_skew("b")

  def test_get_wrong_orientation(self):
    with self.assertRaises(RuntimeError):
      a = self.hist.get_orientation("b")

  def test_does_not_exist(self):
    with self.assertRaises(RuntimeError):
      c = self.hist.get_scalar("g")

class TestActuallyStores(unittest.TestCase):
  def setUp(self):
    self.scalar = 2.5
    self.vector = tensors.Vector(np.array([1.0,2.0,3.0]))
    self.ranktwo = tensors.RankTwo(np.array([
      [1.0,2.0,3.0],
      [4.0,5.0,6.0],
      [7.0,8.0,9.0]]))
    self.symmetric = tensors.Symmetric(np.eye(3))
    q = np.array([[1.0,2.0,3.0],[3.0,4.0,5.0],[6.0,7.0,8.0]])
    q = 0.5 * (q - q.T)
    self.skew = tensors.Skew(q)
    self.orientation = rotations.Orientation(30.0, 60.0, 80.0, angle_type = "degrees")

    self.hist1 = history.History()
    self.hist2 = history.History(store = False)

    self.add_all(self.hist1)
    self.add_all(self.hist2)

    self.storage = np.zeros((self.hist2.size,))
    self.hist2.set_data(self.storage)

    self.set_all(self.hist1)
    self.set_all(self.hist2)

  def add_all(self, hist):
    hist.add_scalar("scalar")
    hist.add_vector("vector")
    hist.add_ranktwo("ranktwo")
    hist.add_symmetric("symmetric")
    hist.add_skew("skew")
    hist.add_orientation("orientation")

  def set_all(self, hist):
    hist.set_scalar("scalar", self.scalar)
    hist.set_vector("vector", self.vector)
    hist.set_ranktwo("ranktwo", self.ranktwo)
    hist.set_symmetric("symmetric", self.symmetric)
    hist.set_skew("skew", self.skew)
    hist.set_orientation("orientation", self.orientation)

  def test_actually_stored(self):
    self.assertTrue(np.isclose(self.scalar, self.storage[0]))

  def test_scalar_stored(self):
    self.assertTrue(np.isclose(self.scalar, self.hist1.get_scalar("scalar")))
    self.assertTrue(np.isclose(self.scalar, self.hist2.get_scalar("scalar")))

  def test_vector_stored(self):
    self.assertEqual(self.vector, self.hist1.get_vector("vector"))
    self.assertEqual(self.vector, self.hist2.get_vector("vector"))

  def test_ranktwo_stored(self):
    self.assertEqual(self.ranktwo, self.hist1.get_ranktwo("ranktwo"))
    self.assertEqual(self.ranktwo, self.hist2.get_ranktwo("ranktwo"))

  def test_symmetric_stored(self):
    self.assertEqual(self.symmetric, self.hist1.get_symmetric("symmetric"))
    self.assertEqual(self.symmetric, self.hist2.get_symmetric("symmetric"))

  def test_skew_stored(self):
    self.assertEqual(self.skew, self.hist1.get_skew("skew"))
    self.assertEqual(self.skew, self.hist2.get_skew("skew"))

  def test_orientation_stored(self):
    self.assertTrue(
        np.allclose(self.orientation.quat, self.hist1.get_orientation("orientation").quat))
    self.assertTrue(
        np.allclose(self.orientation.quat, self.hist2.get_orientation("orientation").quat))

class TestVectorMath(unittest.TestCase):
  def setUp(self):
    self.scalar = 2.5
    self.vector = tensors.Vector(np.array([1.0,2.0,3.0]))
    self.ranktwo = tensors.RankTwo(np.array([
      [1.0,2.0,3.0],
      [4.0,5.0,6.0],
      [7.0,8.0,9.0]]))
    self.symmetric = tensors.Symmetric(np.eye(3))
    q = np.array([[1.0,2.0,3.0],[3.0,4.0,5.0],[6.0,7.0,8.0]])
    q = 0.5 * (q - q.T)
    self.skew = tensors.Skew(q)

    self.hist = history.History()
    self.add_all(self.hist)
    self.set_all(self.hist)

  def add_all(self, hist):
    hist.add_scalar("scalar")
    hist.add_vector("vector")
    hist.add_ranktwo("ranktwo")
    hist.add_symmetric("symmetric")
    hist.add_skew("skew")

  def set_all(self, hist):
    hist.set_scalar("scalar", self.scalar)
    hist.set_vector("vector", self.vector)
    hist.set_ranktwo("ranktwo", self.ranktwo)
    hist.set_symmetric("symmetric", self.symmetric)
    hist.set_skew("skew", self.skew)
  
  def test_add_all(self):
    mult = 2.2
    self.hist.scalar_multiply(mult)
    self.assertTrue(np.isclose(self.hist.get_scalar("scalar"), 
      mult * self.scalar))
    self.assertEqual(self.hist.get_vector("vector"),
        mult * self.vector)
    self.assertEqual(self.hist.get_ranktwo("ranktwo"),
        mult * self.ranktwo)
    self.assertEqual(self.hist.get_symmetric("symmetric"),
        mult * self.symmetric)
    self.assertEqual(self.hist.get_skew("skew"),
        mult * self.skew)

class TestCarefulStore(unittest.TestCase):
  def setUp(self):
    self.hist = history.History(False)
    self.hist.add_scalar("a")
    self.hist.add_scalar("b")
    self.hist.add_scalar("c")
    self.data = np.array([1.0,2.0,3.0])
    self.hist.set_data(self.data)

  def test_got_from_data(self):
    self.assertEqual(self.data[0], self.hist.get_scalar("a"))

  def test_set_data(self):
    self.hist.set_scalar("a", 2.0)
    self.assertEqual(self.hist.get_scalar("a"), 2.0)
    self.assertEqual(self.data[0], 2.0)

class TestSplitStore(unittest.TestCase):
  def setUp(self):
    self.hist = history.History()
    self.hist.add_scalar("a")
    self.hist.set_scalar("a", 1.0)
    self.hist.add_scalar("b")
    self.hist.set_scalar("b", 2.0)
    self.hist.add_scalar("c")
    self.hist.set_scalar("c", 3.0)
    self.hist.add_scalar("d")
    self.hist.set_scalar("d", 4.0)

  def test_reorder(self):
    self.hist.reorder(["b", "c", "a", "d"])
    self.assertEqual(self.hist.items, ["b", "c", "a", "d"])
    self.assertTrue(np.allclose(np.array([2.0,3.0,1.0,4.0]), np.array(self.hist)))

  def test_good_split(self):
    nhist = self.hist.split(["a", "b"])
    self.assertEqual(nhist.size, 2)
    self.assertTrue(np.allclose(np.array([3.0,4.0]), np.array(nhist)))
    self.assertEqual(nhist.items, ["c", "d"])

  def test_bad_split(self):
    with self.assertRaises(RuntimeError):
      nhist = self.hist.split(["b", "c"])

  def test_subset(self):
    subset = self.hist.subset(["a", "d"])
    self.assertTrue(np.allclose(np.array([1.0,4.0]), np.array(subset)))
    self.assertEqual(subset.items, ["a", "d"])

  def test_bad_subset(self):
    with self.assertRaises(RuntimeError):
      subset = self.hist.subset(["a", "no"])

class TestSplitNoStore(unittest.TestCase):
  def setUp(self):
    self.hist = history.History(False)
    self.hist.add_scalar("a")
    self.hist.add_scalar("b")
    self.hist.add_scalar("c")
    self.data = np.array([1.0,2.0,3.0])
    self.hist.set_data(self.data)

  def test_affect(self):
    nhist = self.hist.split(["a"])
    self.assertEqual(nhist.get_scalar("b"), 2.0)
    nhist.set_scalar("b", -2.0)
    self.assertEqual(nhist.get_scalar("b"), -2.0)
    self.assertEqual(nhist.get_scalar("b"), self.hist.get_scalar("b"))
    self.assertEqual(self.data[1], -2.0)

class TestUnion(unittest.TestCase):
  def setUp(self):
    self.scalar1 = 2.0
    self.scalar2 = 3.0
    self.vector1 = tensors.Vector(np.array([1.0,2.0,3.0]))

    self.hist1 = history.History()
    self.hist1.add_scalar("a")
    self.hist1.set_scalar("a", self.scalar1)
    self.hist1.add_vector("b")
    self.hist1.set_vector("b", self.vector1) 

    self.hist2 = history.History()
    self.hist2.add_scalar("c")
    self.hist2.set_scalar("c", self.scalar2)

  def test_union(self):
    hist = self.hist2.add_union(self.hist1)
    self.assertTrue(np.isclose(hist.get_scalar("a"), self.scalar1))
    self.assertTrue(np.isclose(hist.get_scalar("c"), self.scalar2))
    self.assertEqual(hist.get_vector("b"), self.vector1)
