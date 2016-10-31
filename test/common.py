import numpy as np
import itertools

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
