#!/usr/bin/env python3

import sys
sys.path.append('../..')

import numpy as np

from neml.cp import crystallography, slipharden, sliprules, inelasticity, kinematics, singlecrystal
from neml.math import rotations, tensors, nemlmath
from neml import elasticity

import matplotlib.pyplot as plt

if __name__ == "__main__":
  L = np.array([[-0.5,0,0],[0,1.0,0],[0,0,-0.5]])
  erate = 1.0e-4
  steps = 50
  emax = 0.25

  E = 100000.0
  nu = 0.3

  t0_1 = 50.0
  ts_1 = 50.0
  b_1 = 10.0

  t0_2 = 25.0
  ts_2 = -75
  b_2 = 1.0

  g0 = 1.0
  n = 12.0

  # Setup
  L *= erate
  dt = emax / steps / erate

  d = nemlmath.sym(0.5*(L+L.T))
  w = nemlmath.skew(0.5*(L-L.T))

  strengthmodel_1 = slipharden.VoceSlipHardening(ts_1, b_1, t0_1)
  strengthmodel_2 = slipharden.VoceSlipHardening(ts_2, b_2, t0_2)

  strengthmodel = slipharden.SumSlipSingleStrengthHardening([strengthmodel_1,
    strengthmodel_2])

  #strengthmodel = strengthmodel_2

  slipmodel = sliprules.PowerLawSlipRule(strengthmodel, g0, n)
  imodel = inelasticity.AsaroInelasticity(slipmodel)
  emodel = elasticity.IsotropicLinearElasticModel(E, "youngs", nu, "poissons")
  kmodel = kinematics.StandardKinematicModel(emodel, imodel)

  lattice = crystallography.CubicLattice(1.0)
  lattice.add_slip_system([1,1,0],[1,1,1])

  model = singlecrystal.SingleCrystalModel(kmodel, lattice)

  T = 300.0

  h_n = model.init_store()

  s_n = np.zeros((6,))

  d_n = np.zeros((6,))
  w_n = np.zeros((3,))

  u_n = 0.0
  p_n = 0.0
  t_n = 0.0
  T_n = T

  e = [0.0]
  s = [0.0]

  for i in range(steps):
    d_np1 = d_n + d * dt
    w_np1 = w_n + w * dt
    t_np1 = t_n + dt
    T_np1 = T_n

    s_np1, h_np1, A_np1, B_np1, u_np1, p_np1 = model.update_ld_inc(d_np1, d_n, w_np1, w_n, T_np1, T_n, t_np1, t_n, s_n, h_n, u_n, p_n)

    e.append(d_np1[1])
    s.append(s_np1[1])
    
    d_n = np.copy(d_np1)
    w_n = np.copy(w_np1)
    s_n = np.copy(s_np1)
    h_n = np.copy(h_np1)
    t_n = t_np1
    T_n = T_np1
    u_n = u_np1
    p_n = p_np1
  
  plt.plot(e, s)
  plt.show()
