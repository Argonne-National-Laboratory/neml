#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import odf, harmonics, polefigures, crystallography
from neml.math import rotations

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 2000
  O = 5

  orientations = rotations.random_orientations(N)
  
  p = odf.HarmonicsODF(O)
  p.project(orientations)

  lattice = crystallography.CubicLattice(1.0)
  
  for i in range(10):
    print(p.value(rotations.random_orientations(1)[0]))

  polefigures.pole_figure_odf(p, [1,1,1], lattice)
  plt.show()
