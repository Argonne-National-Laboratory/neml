#!/usr/bin/env python3

import sys
sys.path.append('../../..')

from neml.math import rotations
from neml.cp import harmonics

import numpy as np
from scipy.special import legendre, sph_harm

def ds_S2(theta, phi):
  return np.sin(theta)

if __name__ == "__main__":
  print("S2")
  N = 10
  for n in range(N+1):
    theta, phi, wts = harmonics.quadrature_S2(n)
    theta = np.array(theta)
    phi = np.array(phi)
    wts = np.array(wts)
    s = np.sum(wts)
    print("\tOrder %i:" % n)
    print("\t\tVolume ratio %f" % (s / (4.0*np.pi)))

    theta, phi, wts = harmonics.quadrature_S2(2*n)
    theta = np.array(theta)
    phi = np.array(phi)
    wts = np.array(wts)
    s = np.sum(wts)

    for jj in range(-n,n+1):
      for i in range(n+1):
        for j in range(-i,i+1):
          res = np.sum(wts * sph_harm(j,i,theta,phi) * 
              np.conj(sph_harm(jj,n,theta,phi)))
          if (n==i) and (jj==j):
            if not np.isclose(np.real(res),1.0):
              print("\t\t\t(%i,%i)<=>(%i,%i) not 1!" 
                  % (n,jj,i,j))
          else:
            if not np.isclose(np.real(res),0.0):
              print("\t\t\t(%i,%i)<=>(%i,%i) not 0!"
                  % (n,jj,i,j))

  print("SO(3)")
  N = 5
  for n in range(N+1):
    pts, wts = harmonics.quadrature_SO3(n)
    wts = np.array(wts)
    s = np.sum(wts)
    print("\tOrder %i:" % n)
    print("\t\tVolume ratio %f" % (s / (8.0*np.pi**2)))

    pts, wts = harmonics.quadrature_SO3(2*n)
    wts = np.array(wts)
    
    for i in range(-n,n+1):
      for j in range(-n,n+1):
        for nn in range(n+1):
          for ii in range(-nn,nn+1):
            for jj in range(-nn,nn+1):
              val = sum(wt * harmonics.harmonic_SO3(n,i,j,q) * 
                  np.conj(harmonics.harmonic_SO3(nn,ii,jj,q)) for
                    wt,q in zip(wts, pts))
              if ((n==nn) and (i==ii) and (j==jj)):
                if not np.isclose(np.real(val),1,atol=1.0e-4):
                  print("\t\t\t(%i,%i,%i)<=>(%i,%i,%i) not 1 -- %f + %fi!" % (n,i,j,
                      nn,ii,jj,np.real(val),np.imag(val)))
              else:
                if not np.isclose(np.real(val),0,atol=1.0e-4):
                  print("\t\t\t(%i,%i,%i)<=>(%i,%i,%i) not 0 -- %f + %fi!" % (n,i,j,
                    nn,ii,jj,np.real(val),np.imag(val)))

