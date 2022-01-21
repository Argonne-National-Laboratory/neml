#!/usr/bin/env python3

import os.path
import sys
import argparse

import numpy as np

import glob
import unittest

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from regressiontools import rtt, rtt_restart, rtt_pickle

def test_all_regression():
  tests = glob.glob(os.path.join(os.path.dirname(__file__),"test_*"))
  npe = glob.glob(os.path.join(os.path.dirname(__file__),"test*.py"))
  # This removes the files test_*.py, which are unfortunately how I have
  # to name the python files for nose to find them
  for n in npe:
    tests.remove(n)

  for test in tests:
    yield check_regression, test

def test_all_regression_restart():
  tests = glob.glob(os.path.join(os.path.dirname(__file__),"test_*"))
  npe = glob.glob(os.path.join(os.path.dirname(__file__),"test*.py"))
  # This removes the files test_*.py, which are unfortunately how I have
  # to name the python files for nose to find them
  for n in npe:
    tests.remove(n)

  for test in tests:
    yield check_restart_regression, test

def test_all_regression_pickle():
  tests = glob.glob(os.path.join(os.path.dirname(__file__),"test_*"))
  npe = glob.glob(os.path.join(os.path.dirname(__file__),"test*.py"))
  # This removes the files test_*.py, which are unfortunately how I have
  # to name the python files for nose to find them
  for n in npe:
    tests.remove(n)

  for test in tests:
    yield check_restart_pickle, test

def check_regression(tdir):
  print(tdir)
  reference, run = rtt(tdir)

  print(reference[1])
  print(run[1])
  
  assert(np.allclose(reference,run))

def check_restart_regression(tdir):
  print(tdir)
  reference, run = rtt_restart(tdir)

  print(reference[1])
  print(run[1])

  assert(np.allclose(reference, run, atol = 1e-2))

def check_restart_pickle(tdir):
  print(tdir)
  reference, run = rtt_pickle(tdir)

  print(reference[1])
  print(run[1])

  assert(np.allclose(reference, run, atol = 1e-2))
