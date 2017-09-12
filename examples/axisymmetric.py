#!/usr/bin/env python

import sys
sys.path.append('..')

from neml import axisym

import numpy as np

if __name__ == "__main__":
  amodel = axisym.AxisymmetricProblem([10.0,15.0,25.0], [1,1], [10,10], 
      lambda r, t: 0.0, lambda t: 0.0, bias = True)

