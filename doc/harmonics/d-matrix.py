#!/usr/bin/env python3

import numpy as np
import numpy.polynomial as npp
from math import factorial

from sympy import N
from sympy.physics.quantum.spin import Rotation

def reference(l, m, k, x):
  return complex(N(Rotation.d(l, m, k, x).doit()))

if __name__ == "__main__":
  print(reference(0,0,0,np.pi/4))
