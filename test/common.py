import numpy as np
import itertools

mandel = ((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))
mandel_mults = (1,1,1,np.sqrt(2),np.sqrt(2),np.sqrt(2))

def ms2ts(C):
  """
    Convert a Mandel notation stiffness matrix to a full stiffness tensor.
  """
  Ct = np.zeros((3,3,3,3))
  for a in range(6):
    for b in range(6):
      ind_a = itertools.permutations(mandel[a], r=2)
      ind_b = itertools.permutations(mandel[b], r=2)
      ma = mandel_mults[a]
      mb = mandel_mults[b]
      indexes = tuple(ai+bi for ai, bi in itertools.product(ind_a, ind_b))
      for ind in indexes:
        Ct[ind] = C[a,b] / (ma*mb)

  return Ct

def ts2ms(C):
  """
    Convert a stiffness tensor into a Mandel notation stiffness matrix
  """
  Cv = np.zeros((6,6))
  for i in range(6):
    for j in range(6):
      ma = mandel_mults[i]
      mb = mandel_mults[j]
      Cv[i,j] = C[mandel[i]+mandel[j]] * ma * mb

  return Cv

def sym(A):
  """
    Take a symmetric matrix to the Mandel convention vector.
  """
  return np.array([A[0,0], A[1,1], A[2,2], np.sqrt(2)*A[1,2], 
    np.sqrt(2)*A[0,2], np.sqrt(2)*A[0,1]])

def usym(v):
  """
    Take a Mandel symmetric vector to the full matrix.
  """
  return np.array([
    [v[0], v[5]/np.sqrt(2), v[4]/np.sqrt(2)],
    [v[5]/np.sqrt(2), v[1], v[3]/np.sqrt(2)],
    [v[4]/np.sqrt(2), v[3]/np.sqrt(2), v[2]]
    ])

def skew(A):
  """
    Take a skew matrix to my vector convention.
  """
  return np.array([-A[1,2], A[0,2], -A[0,1]])

def uskew(v):
  """
    Take a skew vector in my convention to the matrix
  """
  return np.array([
    [0, -v[2], v[1]],
    [v[2], 0, -v[0]],
    [-v[1], v[0], 0]
    ])

def differentiate(fn, x0, eps = 1.5e-6):
  if np.isscalar(x0):
    sx = (1,)
  else:
    sx = x0.shape

  f0 = fn(x0)

  if np.isscalar(f0):
    sf = (1,)
  else:
    sf = f0.shape

  tshape = sf + sx
  rshape = sx + sf

  D = np.zeros(rshape)

  for index in np.ndindex(sx):
    if np.isscalar(x0):
      diff = np.abs(x0 * eps)
    else:
      diff = np.abs(x0[index] * eps)
    if diff < eps:
      diff = eps
    
    if np.isscalar(x0):
      xp = diff
    else:
      xp = np.zeros(sx)
      xp[index] = diff

    fp = fn(xp+x0)

    D[index] = (fp - f0) / diff
  
  Df = np.zeros(tshape)
  # Reverse
  for ind1 in np.ndindex(sx):
    for ind2 in np.ndindex(sf):
      Df[ind2+ind1] = D[ind1+ind2]
  
  if np.isscalar(x0):
    return Df[0][0]
  else:
    return Df

def make_dev(s):
  return s - np.array([1,1,1,0,0,0]) * np.sum(s[:3]) / 3.0
