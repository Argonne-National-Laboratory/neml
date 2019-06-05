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

class TestCopy(unittest.TestCase):
  def setUp(self):
    self.scalar = 2.5

    self.hist = history.History()
    self.hist.add_scalar("a")
    self.storage = np.zeros((self.hist.size,))
    self.hist.set_scalar("a", self.scalar)

  def should_change(self):
    self.storage[0] = 3.0
    self.assertTrue(np.isclose(self.hist.get_scalar("a"), 3.0))

  def should_not_change(self):
    newhist = self.hist.deepcopy()
    self.storage[0] = 3.0
    self.assertFalse(np.isclose(newhist.get_scalar("a"), self.hist.get_scalar("a")))

    self.hist.set_scalar("a", 4.0)
    self.assertFalse(np.isclose(newhist.get_scalar("a"), self.hist.get_scalar("a")))

    self.assertTrue(np.isclose(self.hist.get_scalar("a"), 4.0))
    self.assertTrue(np.isclose(newhist.get_scalar("a"), self.scalar))
