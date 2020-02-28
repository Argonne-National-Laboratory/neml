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
    pass

class TestMatDefinition(unittest.TestCase):
  def setUp(self):
    pass
