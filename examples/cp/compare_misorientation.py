#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

import time

from neml.cp import crystallography
from neml.math import rotations

def loop(A, B, sgroup):
  return [sgroup.misorientation(a,b) for a,b in zip(A,B)]

if __name__ == "__main__":
  N = 50000

  ref1 = rotations.Orientation(30.0,60.0,45.0, angle_type = "degrees")
  ref2 = rotations.Orientation(1.2,55.0,33.0, angle_type = "degrees")

  A = rotations.random_orientations(N)
  B = rotations.random_orientations(N)

  sgroup = crystallography.SymmetryGroup("432")

  print(sgroup.misorientation(ref1, ref2))
  print(sgroup.misorientation_block([ref1], [ref2])[0])

  ta = time.time()
  res1 = loop(A, B, sgroup)
  tb = time.time()

  t1 = time.time()
  res2 = sgroup.misorientation_block(A, B)
  t2 = time.time()

  print("Loop method: %f seconds" % (tb-ta))
  print("Block method: %f seconds" % (t2-t1))
