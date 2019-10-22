from neml.math import tensors, rotations
from neml.cp import crystallography

import numpy as np
import numpy.linalg as la

from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

limit_lambert_equal_area = np.sqrt(2)
limit_stereographic = 1.0

triangle_limits = {"cubic": (0,0)}

def project_stereographic(v):
  """
    Stereographic projection of the given vector into a numpy array

    Parameters:
      v:        input vector
  """
  return np.array([v[0]/(1.0+v[2]), v[1]/(1.0+v[2])])

def project_lambert_equal_area(v):
  """
    Lambert equal area projection of the given vector into a numpy array

    Parameters:
      v:        input vector
  """
  return np.array([np.sqrt(2/(1+v[2]))*v[0]/np.sqrt(2),
    np.sqrt(2/(1+v[2]))*v[0]/np.sqrt(2)])

def inverse_project_stereographic(pt):
  """
    Inverse stereographic projection

    Parameters:
      pt:            point, as numpy array
  """
  X = pt[0]
  Y = pt[1]
  bot = 1.0 + X**2.0 + Y**2.0
  return tensors.Vector([2.0*X/bot,2.0*Y/bot,(X**2.0+Y**2.0-1.0)/bot])

def inverse_project_lambert_equal_area(pt):
  """
    Inverse Lambert projections

    Parameters:
      pt:            point, as a numpy array
  """
  X = pt[0]
  Y = pt[1]
  f = np.sqrt(1.0-(X**2.0+Y**2.0)/4)
  return tensors.Vector([f*X,f*Y,-1.0+(X**2.0+Y**2.0)/2])

def cart2pol(x):
  """
    Convert a cartesian point into polar coordinates

    Parameters:
      x:        point to convert to polar
  """
  return np.array([np.arctan2(x[1], x[0]), la.norm(x)])

def pole_figure_discrete(orientations, pole, lattice,  
    projection = "stereographic", 
    sample_symmetry = crystallography.symmetry_rotations("222"),
    x = tensors.Vector([1.0,0,0]), y = tensors.Vector([0,1.0,0]),
    axis_labels = ["X", "Y"]):
  """
    Plot a standard pole figure given a collection of discrete points.

    Parameters:
      orientations:      list of orientations
      pole:              pole as a integer list
      lattice:           crystal lattice class, including crystal symmetry

    Keyword Args:
      projection:        which projection to use
      sample_symmetry:   what sample symmetry to apply to the pole figure,
                         defaults to orthorhombic
      x:                 x direction for plot, defaults to [1,0,0]
      y:                 y direction for plot, defaults to [0,1,0]
      axis_labels:       axis labels to include on the figure
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
    lim = limit_stereographic
  elif projection == "lambert":
    pop = project_lambert_equal_area
    lim = limit_lambert_equal_area
  else:
    raise ValueError("Unknown projection %s" % projection)

  cart = np.array([pop(v) for v in pts])

  # Convert to polar
  polar = np.array([cart2pol(cp) for cp in cart])

  # Plot
  ax = plt.subplot(111, projection='polar')
  ax.scatter(polar[:,0], polar[:,1], c='k', s=10.0)

  # Make the graph nice
  plt.ylim([0,lim])
  ax.grid(False)
  ax.get_yaxis().set_visible(False)
  if len(axis_labels) == 2:
    plt.xticks([0,np.pi/2], axis_labels)
  else:
    ax.get_xaxis().set_visible(False)
  ax.xaxis.set_minor_locator(plt.NullLocator())

def pol2cart(R, T):
  """
    Convert polar coordinates to cartesian

    Parameters:
      R:         radii
      T:         angle, in radians
  """
  X = R * np.cos(T)
  Y = R * np.sin(T)

  return X, Y

def inverse_pole_figure_discrete(orientations, direction, lattice,  
    reduce_figure = False, color = False,
    sample_symmetry = crystallography.symmetry_rotations("222"),
    x = [1,0,0], y = [0,1,0], axis_labels = None, nline = 100):
  """
    Plot an inverse pole figure given a collection of discrete points.

    Parameters:
      orientations:      list of orientations
      direction:         sample direction
      lattice:           crystal lattice class, including crystal symmetry

    Keyword Args:
      reduce_figure:     reduce to some fundamental region.
                         Options include
                         False -- don't do it
                         "cubic" -- cubic convention
                         [v1,v2,v3] -- list of crystallographic points defining the triangle
      color:             color points based on the provided triangle
      sample_symmetry:   what sample symmetry to apply to the pole figure,
                         defaults to orthorhombic
      x:                 crystallographic x direction for plot, defaults to [1,0,0]
      y:                 crystallographic y direction for plot, defaults to [0,1,0]
      axis_labels:       axis labels to include on the figure
      nline:             number of discrete points to use in plotting lines on the triangle
  """
  pts = np.vstack(tuple(project_ipf(q, lattice, direction, 
    sample_symmetry = sample_symmetry, x = x, y = y) for q in orientations))

  if reduce_figure:
    if reduce_figure == "cubic":
      vs = (np.array([0,0,1.0]), np.array([1.0,0,1]), np.array([1.0,1,1]))
    elif len(reduce_figure) == 3:
      vs = reduce_figure
    else:
      raise ValueError("Unknown reduction type %s!" % reduce_figure)
    
    pts = reduce_points_triangle(pts, v0 = vs[0], v1=vs[1], v2=vs[2])

  pop = project_stereographic
  lim = limit_stereographic

  cpoints = np.array([pop(v) for v in pts])

  # Make the graph nice
  if reduce_figure:
    ax = plt.subplot(111)
    if color:
      rgb = ipf_color(pts, v0 = vs[0], v1 = vs[1], v2=vs[2]) 
      ax.scatter(cpoints[:,0], cpoints[:,1], c=rgb, s = 10.0)
    else:
      ax.scatter(cpoints[:,0], cpoints[:,1], c='k', s = 10.0) 
    ax.axis('off')
    if axis_labels:
      plt.text(0.1,0.11,axis_labels[0], transform = plt.gcf().transFigure)
      plt.text(0.86,0.11,axis_labels[1], transform = plt.gcf().transFigure)
      plt.text(0.74,0.88,axis_labels[2], transform = plt.gcf().transFigure)
    
    for i,j in ((0,1),(1,2),(2,0)):
      v1 = vs[i]
      v2 = vs[j]
      fs = np.linspace(0,1,nline)
      pts = np.array([pop((f*v1+(1-f)*v2)/la.norm(f*v1+(1-f)*v2)) for f in fs])
      plt.plot(pts[:,0], pts[:,1], color = 'k')
  else:
    polar = np.array([cart2pol(v) for v in cpoints])
    ax = plt.subplot(111, projection='polar')
    ax.scatter(polar[:,0], polar[:,1], c='k', s=10.0)
    plt.ylim([0,lim])
    ax.grid(False)
    ax.get_xaxis().set_visible(False)    
    ax.get_yaxis().set_visible(False)
    ax.xaxis.set_minor_locator(plt.NullLocator())

def project_ipf(q, lattice, direction, 
    sample_symmetry = crystallography.symmetry_rotations("222"), 
    x = [1,0,0], y = [0,1,0]):
  """
    Project a single sample direction onto a crystal

    Parameters:
      q:                 lattice orientation
      lattice:           lattice object describing the crystal system
      direction:         sample direction 
    
    Keyword Args:
      sample_symmetry:   sample symmetry operators
      x:                 x direction (crystallographic) of the projection
      y:                 y direction (crystallographic) of the projection
  """
  xv = lattice.miller2cart_direction(x).normalize()
  yv = lattice.miller2cart_direction(y).normalize()

  if not np.isclose(xv.dot(yv),0.0):
    raise ValueError("Lattice directions are not orthogonal!")

  zv = xv.cross(yv)
  
  trans = rotations.Orientation(np.vstack((xv.data,yv.data,zv.data)))
  
  d = tensors.Vector(direction).normalize()

  pts = []
  for srot in sample_symmetry:
    # Sample directions
    spt = srot.apply(d)
    
    # Crystal coordinates
    cpt = q.apply(spt)

    # Lattice symmetry
    for op in lattice.symmetry.ops:
      cppt = op.apply(cpt)

      # Into the right coordinates
      fpt = trans.apply(cppt)

      pts.append(fpt.data)
  
  # By convention just keep the points in the upper hemisphere
  pts = np.array(pts)
  pts = pts[pts[:,2]>0]
  
  return pts

def reduce_points_triangle(pts, v0 = np.array([0,0,1.0]), 
    v1 = np.array([1.0,0,1]), v2 = np.array([1.0,1,1])):
  """
    Reduce points to a standard stereographic triangle

    Default vectors give you the cubic triangle

    Parameters:
      pts:       projected points

    Keyword Args:
      v1:        1st vector
      v2:        2nd vector
      v3:        3rd vector
  """
  v0 /= la.norm(v0)
  v1 /= la.norm(v1)
  v2 /= la.norm(v2)

  n0 = np.cross(v0, v1)
  n1 = np.cross(v1, v2)
  n2 = np.cross(v2, v0)

  npts = []
  for pt in pts:
    if np.dot(pt, n0) >= 0 and np.dot(pt, n1) >= 0 and np.dot(pt, n2) >= 0:
      npts.append(pt)

  return np.array(npts)

def ipf_color(pts, v0 = np.array([0,0,1.0]), 
    v1 = np.array([1.0,0,1]), v2 = np.array([1.0,1,1])):
  """
    Color ipf points based on the provided triangle

    Assumes the points already fall inside the appropriate triangle

    Parameters:
      pts:       list of points

    Keyword Args:
      v0:        point on triangle
      v1:        point on triangle
      v2:        point on triangle
  """
  v0 /= la.norm(v0)
  v1 /= la.norm(v1)
  v2 /= la.norm(v2)

  vs = [v0, v1, v2]

  colors = np.zeros((len(pts),3))

  for i,pt in enumerate(pts):
    for j in range(3):
      colors[i,j] = np.dot(vs[j], pt)

  for i in range(3):
    mf = 1.0
    for j in range(3):
      t = np.dot(vs[i], vs[j])
      if t < mf:
        mf = t
    
    colors[:,i] -= mf
    colors[:,i] /= (1.0 - mf)

  return colors

def ipf_color_chart(v0 = np.array([0,0,1.0]), 
    v1 = np.array([1.0,0,1]), v2 = np.array([1.0,1,1]), 
    axis_labels = None, ngrid = 50, nline = 100):
  """
    Make a color chart for the IPF coloring used here

    Keyword Args:
      v0:            point on triangle
      v1:            point on triangle
      v2:            point on triangle
      axis_labels:   labels to draw on triangle
      ngrid:         number of grid points to use
      nline:         number of line points to use
  """
  vs = [v0,v1,v2]

  pop = project_stereographic
  lim = limit_stereographic
  def ipop(x):
    x1 = inverse_project_stereographic(x).data
    if x1[2] < 0:
      x1[2] *= -1
    return x1

  ax = plt.subplot(111)
  ax.axis('off')
  if axis_labels:
    plt.text(0.12,0.11,axis_labels[0], transform = plt.gcf().transFigure)
    plt.text(0.86,0.11,axis_labels[1], transform = plt.gcf().transFigure)
    plt.text(0.74,0.90,axis_labels[2], transform = plt.gcf().transFigure)

  pts = np.array([pop(v0/la.norm(v0)), pop(v1/la.norm(v1)), pop(v2/la.norm(v2))])
  xrang = np.linspace(np.min(pts[:,0]), np.max(pts[:,0]), ngrid)
  yrang = np.linspace(np.min(pts[:,1]), np.max(pts[:,1]), ngrid)

  X, Y = np.meshgrid(xrang,yrang)
  colors = np.zeros(X.shape+(3,))
  for ind in np.ndindex(X.shape):
    colors[ind] = ipf_color([ipop([X[ind],Y[ind]])], v0 = v0, v1 = v1, v2 = v2)[0]
  
  im = plt.imshow(colors, extent = [xrang[0], xrang[-1], yrang[0], yrang[-1]], 
      origin = 'lower')
  
  net_pts = []
  for i,j in ((0,1),(1,2),(2,0)):
    vi = vs[i]
    vj = vs[j]
    fs = np.linspace(0,1,nline)
    pts = np.array([pop(((1-f)*vi+(f)*vj)/la.norm((1-f)*vi+(f)*vj)) for f in fs])
    net_pts.extend(list(pts[:-1]))

  net_pts = np.array(net_pts)
  pat = Polygon(net_pts, closed = True, edgecolor='k', fc = None, fill = False) 
  ax.add_artist(pat)
  im.set_clip_path(pat)
