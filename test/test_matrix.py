#!/usr/bin/env python

from neml.math import matrix

import common

import unittest
import numpy as np

class TestVector(unittest.TestCase):
  def setUp(self):
    self.a_data = [1.0,1.1,2.3]
    self.a = matrix.FlatVector(self.a_data)

  def test_definition(self):
    self.assertTrue(np.allclose(self.a_data, self.a))

  def test_size(self):
    self.assertEqual(self.a.n, len(self.a_data))

class TestBasicMatrix(unittest.TestCase):
  def setUp(self):
    self.np_mat = np.array([[2.0,3.0,1.0],[7.1,3.1,-1.0],[2.0,-1.2,8.1]])
    self.M = matrix.SquareMatrix(3, type = "dense", 
        data = list(self.np_mat.flatten())) 

  def test_definition(self):
    self.assertTrue(np.allclose(self.np_mat, self.M))

  def test_size(self):
    self.assertEqual(self.M.m, self.np_mat.shape[0])
    self.assertEqual(self.M.n, self.np_mat.shape[1])

class TestMatVec(unittest.TestCase):
  def setUp(self):
    self.a_data = [1.0,1.1,2.3]
    self.a = matrix.FlatVector(self.a_data)

    self.b_data = [-8.0,1.1,2.2,3.3]
    self.b = matrix.FlatVector(self.b_data)

    self.M_data = np.array([[2.0,3.0,1.0],[7.1,3.1,-1.0],[2.0,-1.2,8.1]])
    self.M = matrix.SquareMatrix(3, type = "dense", 
        data = list(self.M_data.flatten())) 

  def test_correct(self):
    self.assertTrue(np.allclose(self.M.dot(self.a), np.dot(self.M_data, 
      self.a_data)))

  def test_bad_size(self):
    with self.assertRaises(ValueError):
      self.M.dot(self.b)

class TestMatDefinition(unittest.TestCase):
  def test_zero(self):
    self.assertTrue(np.allclose(matrix.SquareMatrix(4), np.zeros((4,4))))

  def test_identity(self):
    self.assertTrue(np.allclose(matrix.SquareMatrix(4, type = "identity"), 
      np.eye(4)))

  def test_diagonal(self):
    self.assertTrue(np.allclose(matrix.SquareMatrix(4, type = "diagonal",
      data = [1,2,3,4]), np.diag([1,2,3,4])))

  def test_block_diagonal(self):
    self.assertTrue(np.allclose(matrix.SquareMatrix(4, type = "diagonal_blocks",
      data = [1,2], blocks = [3,1]), np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,2]])))

  def test_blocks(self):
    self.assertTrue(np.allclose(matrix.SquareMatrix(4, type = "block",
      data = [1,2,3,4], blocks=[3,1]),
      np.array([
        [1,1,1,2],
        [1,1,1,2],
        [1,1,1,2],
        [3,3,3,4]])))
