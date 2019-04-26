#!/usr/bin/env python

from neml.math import tensors

import unittest
import numpy as np
import numpy.linalg as la

class TestVector(unittest.TestCase):
  def setUp(self):
    self.a = np.array([2.2,-1.2,2.5])
    self.b = np.array([0.0,5.8,1.1])

    self.va = tensors.Vector(self.a)
    self.vb = tensors.Vector(self.b)

    self.s = 2.1
  
  def test_assign(self):
    self.vb = self.va
    self.assertTrue(np.allclose(self.vb.data, self.a))

  def test_norm(self):
    self.assertAlmostEqual(self.va.norm(), la.norm(self.a))

  def test_dot(self):
    self.assertAlmostEqual(np.dot(self.a,self.b), self.va.dot(self.vb))

  def test_smultiply(self):
    self.assertTrue(np.allclose(
      (self.s * self.va).data,
      self.s * self.a))
    self.assertTrue(np.allclose(
      (self.va * self.s).data,
      self.s * self.a))

  def test_sdivide(self):
    self.assertTrue(np.allclose(
      (self.va / self.s).data,
      self.a / self.s))

  def test_add(self):
    self.assertTrue(np.allclose(
      (self.va + self.vb).data,
      self.a + self.b))
    self.assertTrue(np.allclose(
      (self.vb + self.va).data,
      self.a + self.b))

  def test_subtract(self):
    self.assertTrue(np.allclose(
      (self.va - self.vb).data,
      self.a - self.b))
    self.assertTrue(np.allclose(
      (self.vb - self.va).data,
      self.b - self.a))

  def test_negate(self):
    self.assertTrue(np.allclose(
      (-self.va).data,
      -self.a))
    self.assertTrue(np.allclose(
      self.va.opposite().data,
      -self.a))

  def test_normalize(self):
    self.va.normalize()
    self.assertTrue(np.allclose(
      self.va.data, self.a / la.norm(self.a)))

  def test_cross(self):
    self.assertTrue(np.allclose(
      self.va.cross(self.vb).data,
      np.cross(self.a, self.b)))

  def test_equality(self):
    self.assertTrue(self.va == self.va)
    self.assertFalse(self.va == self.vb)

  def test_inequality(self):
    self.assertFalse(self.va != self.va)
    self.assertTrue(self.va != self.vb)
