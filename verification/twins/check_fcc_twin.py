#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np
import matplotlib.pyplot as plt
from neml.math import rotations
from neml.cp import crystallography, polefigures

if __name__ == "__main__":
  a = 1.3
  lattice = crystallography.CubicLattice(a)  
  lattice.add_twin_system([1,1,2],[1,1,1],[1,1,2],[1,1,-1])

  orientations = [lattice.reorientation(0,i) for i in range(lattice.nslip(0))]
  I = rotations.Orientation(np.array([1,0,0,0.0]))

  for o in orientations:
    m = lattice.symmetry.misorientation(o, I)
    axis, angle = m.to_axis_angle("degrees")
    print(angle)

  polefigures.pole_figure_discrete(orientations,[1,1,1], lattice)
  plt.show()
