#!/usr/bin/env python3

import sys
sys.path.append('../../..')

import numpy as np

from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal, polycrystal, crystaldamage
from neml.math import rotations, tensors, nemlmath
from neml import elasticity, drivers

import matplotlib.pyplot as plt

nthreads = 4

def drive(model, orientations, erate, steps, L, T = 300):
  pmodel = polycrystal.TaylorModel(model, orientations, nthreads = nthreads)

  L *= erate
  dt = emax / steps / erate

  h_n = pmodel.init_store()
 
  d_inc = nemlmath.sym(0.5*(L+L.T))
  w_inc = nemlmath.skew(0.5*(L-L.T))
  
  s_n = np.zeros((6,))
  d_n = np.zeros((6,))
  w_n = np.zeros((3,))

  u_n = 0.0
  p_n = 0.0

  D = [np.zeros(6,)]
  S = [np.zeros(6,)]

  for i in range(steps):
    print(i)
    d_np1 = d_n + d_inc * dt
    w_np1 = w_n + w_inc * dt

    try:
      s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = pmodel.update_ld_inc(
          d_np1, d_n, w_np1, w_n, T, T, dt, 0,
          s_n, h_n, u_n, p_n)
    except RuntimeError:
      break
    
    d_n = np.copy(d_np1)
    w_n = np.copy(w_np1)

    s_n = np.copy(s_np1)
    h_n = np.copy(h_np1)

    u_n = u_np1
    p_n = p_np1

    D.append(d_np1)
    S.append(s_np1)

  return np.array(D), np.array(S)

if __name__ == "__main__":
  emax = 0.5
  N = 25
  erate = 1.0
  steps = 500

  orientations = rotations.random_orientations(N)

  L = np.array([[-0.5,0,0],[0,-0.5,0],[0,0,1.0]])

  lattice = crystallography.CubicLattice(1.0)
  lattice.add_slip_system([1,1,0],[1,1,1])

  nslip = lattice.ntotal

  s0 = 40.0
  k = 100.0
  sat = 100.0
  m = 1.0

  g0 = 1.0
  n = 6.0

  mu = 29000.0
  E = 120000.0
  nu = 0.3

  c = 40
  beta = 3.0

  emodel = elasticity.CubicLinearElasticModel(E, nu, mu, "moduli")

  strengthmodel = slipharden.VocePerSystemHardening(
      [s0] * nslip, [k] * nslip, [sat] * nslip, [m] * nslip)

  slipmodel = sliprules.PowerLawSlipRule(strengthmodel, g0, n)

  imodel = inelasticity.AsaroInelasticity(slipmodel)

  base_kin = kinematics.StandardKinematicModel(emodel, imodel)
  base_model = singlecrystal.SingleCrystalModel(base_kin, lattice)
   
  dmodel = crystaldamage.WorkPlaneDamage()
  func = crystaldamage.SigmoidTransformation(c, beta)
  dmg_odel = crystaldamage.PlanarDamageModel(dmodel, func, func, lattice)
  dmg_kin = kinematics.DamagedStandardKinematicModel(emodel, imodel, dmg_odel)
  dmg_model = singlecrystal.SingleCrystalModel(dmg_kin, lattice)

  D1, S1 = drive(base_model, orientations, erate, steps, L)
  D2, S2 = drive(dmg_model, orientations, erate, steps, L)

  plt.plot(D1[:,2], S1[:,2])
  plt.plot(D2[:,2], S2[:,2])
  plt.show()
