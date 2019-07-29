from neml.math import tensors, rotations
from neml.cp import crystallography

import numpy as np
import numpy.linalg as la

import matplotlib.pyplot as plt

def project_stereographic(v):
  """
    Stereographic projection of the given vector into a numpy array
  """
  return np.array([v[0]/(1.0+v[2]), v[1]/(1.0+v[2])])

def project_lambert_equal_area(v):
  """
    Lambert equal area projection of the given vector into a numpy array
  """
  return np.array([np.sqrt(2/(1+v[2]))*v[0]/np.sqrt(2),
    np.sqrt(2/(1+v[2]))*v[0]/np.sqrt(2)])

def cart2pol(x):
  """
    Convert a cartesian point into polar coordinates
  """
  return np.array([np.arctan2(x[1], x[0]), la.norm(x)])

def pole_figure_discrete(orientations, pole, lattice,  
    projection = "stereographic", 
    sample_symmetry = crystallography.symmetry_rotations("222"),
    x = tensors.Vector([1.0,0,0]), y = tensors.Vector([0,1.0,0]),
    axis_labels = ["X", "Y"]):
  """
    Plot a standard pole figure.

    Parameters:
      orientations      list of orientations
      pole              pole as a integer list
      lattice           crystal lattice class, including crystal symmetry

    Optional:
      projection        which projection to use
      sample_symmetry   what sample symmetry to apply to the pole figure,
                        defaults to orthorhombic
      x                 x direction for plot, defaults to [1,0,0]
      y                 y direction for plot, defaults to [0,1,0]
      axis_labels       axis labels to include on the figure
  """
  # Get the rotation from the standard coordinates to the given x and y
  srot = rotations.Orientation(x,y)
  
  # Get all equivalent poles
  p = lattice.miller2cart_direction(pole)
  eq_poles = lattice.equivalent_vectors(p)
  eq_poles = [pp.normalize() for pp in eq_poles]

  # Apply sample symmetric
  eq_poles = [op.apply(p) for p in eq_poles for op in sample_symmetry]
  
  # Get the points on the sphere
  pts = [srot.apply(o.inverse().apply(pp)) for pp in eq_poles for o in orientations]

  # Get rid of those that are in the lower hemisphere
  pts = [p for p in pts if p[2] > 0.0]

  # Project
  if projection == "stereographic":
    pop = project_stereographic
  elif projection == "lambert":
    pop = project_lambert_equal_area
  else:
    raise ValueError("Unknown projection %s" % projection)

  cart = np.array([pop(v) for v in pts])

  # Convert to polar
  polar = np.array([cart2pol(cp) for cp in cart])

  # Plot
  ax = plt.subplot(111, projection='polar')
  ax.scatter(polar[:,0], polar[:,1], c='k', s=10.0)

  # Make the graph nice
  plt.ylim([0,1])
  ax.grid(False)
  ax.get_yaxis().set_visible(False)
  if len(axis_labels) == 2:
    plt.xticks([0,np.pi/2], axis_labels)
  else:
    ax.get_xaxis().set_visible(False)
  ax.xaxis.set_minor_locator(plt.NullLocator())
