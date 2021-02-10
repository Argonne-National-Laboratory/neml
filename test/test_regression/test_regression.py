#!/usr/bin/env python3

import os.path
import sys
import argparse

import numpy as np

import glob
import unittest

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from regressiontools import rtt

def test_all_regression():
  tests = glob.glob(os.path.join(os.path.dirname(__file__),"test*"))
  npe = glob.glob(os.path.join(os.path.dirname(__file__),"test*.py"))
  for n in npe:
    tests.remove(n)

  for test in tests:
    yield check_regression, test

def check_regression(tdir):
  print(tdir)
  reference, run = rtt(tdir)

  print(reference[1])
  print(run[1])
  
  assert(np.allclose(reference,run))
