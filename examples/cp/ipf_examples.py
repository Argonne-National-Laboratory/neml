#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import polefigures, crystallography
from neml.math import rotations

import matplotlib.pyplot as plt

if __name__ == "__main__":
  N = 1000

  lattice = crystallography.CubicLattice(1.0)

  orientations = rotations.random_orientations(N)
  d = np.array([0,0,1.0])
  
  plt.figure()
  polefigures.inverse_pole_figure_discrete(orientations, d, lattice)
  plt.show()

  plt.figure()
  polefigures.inverse_pole_figure_discrete(orientations, d, lattice, 
      reduce_figure = "cubic", axis_labels = ["100","110","111"])
  plt.show()

  plt.figure()
  polefigures.inverse_pole_figure_discrete(orientations, d, lattice, 
      reduce_figure = "cubic", color = True, axis_labels = ["100","110","111"])
  plt.show()

  plt.figure()
  polefigures.ipf_color_chart(axis_labels = ["100", "110", "111"])
  plt.show()
