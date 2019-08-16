#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import odf, harmonics
from neml.math import rotations

if __name__ == "__main__":
  N = 2000
  O = 5

  orientations = rotations.random_orientations(N)
  
  p = odf.HarmonicsODF(O)
  p.project(orientations)
  
  for i in range(10):
    print(p.value(rotations.random_orientations(1)[0]))
